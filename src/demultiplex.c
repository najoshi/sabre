
/*
 * set softtab=4
 * set shiftwidth=4
 * set expandtab
 *
 */

/* 
 * sabre FASTQ files demultiplexing
 * demultiplex.c: FASTQ demultiplexing
 * 
 */

#include "demultiplex.h"

void* demult_runner(void *arg)
{

    kseq_t *fqrec1;
    kseq_t *fqrec2;
    barcode_data_t *curr;

    int l1, l2;
    thread_data_t* thread_data = (thread_data_t*)arg;
    int my_line_num;

    fqrec1 = kseq_init(thread_data->params->fq1_fd);

    if(thread_data->params->paired > 0) {
        fqrec2 = kseq_init(thread_data->params->fq2_fd);
    }

    /* Get reads, one at a time */

    while(1) {

        // lock reading
        pthread_mutex_lock(thread_data->in_lock);

        l1 = kseq_read(fqrec1);

	// sanity check no more reads
	if(l1 < 0 ) {
            pthread_mutex_unlock(thread_data->in_lock);
            break;	
	}
	
        int n_crop = 0;

        char *actl_bc = NULL;

        char *fqread1 = NULL;
        char *fqread2 = NULL;

        size_t fq_size = 0;

        fq_size += strlen(fqrec1->seq.s);
        fq_size += (strlen(fqrec1->name.s)*2);
        fq_size += strlen(fqrec1->qual.s);
        fq_size += (strlen(fqrec1->comment.s)*2);
        fq_size += 2;// header signs @ and +
        fq_size += 2;//two colons (:)
        fq_size += 4;//cariage returns
        fq_size += 2;//two spaces
        fq_size += 1000;//test

        if(thread_data->params->paired > 0 || thread_data->params->combine > 0) {
            l2 = kseq_read(fqrec2);
            if (l2 < 0) {
                fprintf (stderr, "\n\
                                  \n ERROR: R2 file is shorter than R1 file.\
                                  \n Stopping here:\
                                  \n %s\
                                  \n",
                                  fqrec1->name.s);
                break;
            }
            fq_size += strlen(fqrec2->seq.s);
        }

        /* Step 1: Find matching barcode */
        curr = thread_data->curr;
        while(curr) {
            n_crop = chk_bc_mtch(curr->bc, fqrec1->seq.s, thread_data->params->mismatch, thread_data->params->max_5prime_crop);
            if(n_crop >= 0) {
                //found matching barcode
                actl_bc = strndup( (fqrec1->seq.s)+n_crop, strlen(curr->bc) );
                break;
            }
            curr = thread_data->curr->next;
        }

        // unlock reading
        my_line_num = *(thread_data->line_num);
        *thread_data->line_num += 1;
        pthread_mutex_unlock(thread_data->in_lock);

        /* Step 2: Write read out into barcode specific file */

        // lock writing
        while(*(thread_data->out_line_num) != my_line_num) {
            pthread_cond_wait(thread_data->cv, thread_data->out_lock);
        }
        *thread_data->out_line_num += 1;

        pthread_cond_broadcast(thread_data->cv);  // Tell everyone it might be their turn!

        char *umi_idx = NULL;

        if(curr != NULL) {
            //for now assume barcode and umi are in R1 read
            if(thread_data->params->umi > 0) {

                const char *actl_umi_idx = (fqrec1->seq.s)+strlen(curr->bc)+n_crop;

                if(strlen(actl_umi_idx) < thread_data->params->min_umi_len) {
			//protect by mutex umis_2_short_file
                    fprintf(thread_data->params->umis_2_short_fd, "%s\t%s\t%zu\t%d\n", fqrec1->name.s, actl_umi_idx, strlen(actl_umi_idx), thread_data->params->min_umi_len);
                    continue;
                }
                else {
                   umi_idx = strdup(actl_umi_idx);
                   umi_idx[thread_data->params->min_umi_len] = '\0';
                   fq_size += strlen(umi_idx);
                }
            }

            if(thread_data->params->combine > 0) {
                fqread1 = (char*) malloc(fq_size);
                fqread1[0] = '\0';

                get_merged_fqread(&fqread1, fqrec1, fqrec2, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);
			//protect by mutex umis_2_short_file
                gzwrite(curr->bcfile1, fqread1, strlen(fqread1));
            }
            else {
                fqread1 = (char*) malloc(fq_size + 1);
                fqread2 = (char*) malloc(fq_size + 1);

                fqread1[0] = '\0';
                fqread2[0] = '\0';

                get_fqread(&fqread1, fqrec1, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);
                gzwrite(curr->bcfile1, fqread1, strlen(fqread1));

                if(thread_data->params->paired > 0) {
                    get_fqread(&fqread2, fqrec1, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);
                    //fprintf(curr->bcfile2, "%s", fqread2);
                    gzwrite(curr->bcfile2, fqread2, strlen(fqread2));
                    curr->num_records += 1;
                }
            }
            curr->num_records += 1;
        }
        else {
            fqread1 = (char*) malloc(fq_size + 1);
            fqread2 = (char*) malloc(fq_size + 1);

            fqread1[0] = '\0';
            fqread2[0] = '\0';

            get_fqread(&fqread1, fqrec1, NULL, NULL, thread_data->params->no_comment, 0);
            gzwrite(thread_data->params->unassigned1_fd, fqread1, strlen(fqread1));
            thread_data->metrics->num_unknown += 1;

            if(thread_data->params->paired > 0) {
                get_fqread(&fqread2, fqrec2, NULL, NULL, thread_data->params->no_comment, 0);
                gzwrite(thread_data->params->unassigned2_fd, fqread2, strlen(fqread2));
                thread_data->metrics->num_unknown += 1;
            }
        }

        thread_data->metrics->total += 2;

        // unlock writing
        pthread_mutex_unlock(thread_data->out_lock);

        free(fqread1);
        free(fqread2);
        free(actl_bc);
        free(umi_idx);
    }

    free(thread_data);

    return NULL;
}
