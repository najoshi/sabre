
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

    fq_rec_t *fq_rec1;
    fq_rec_t *fq_rec2;

    fq_rec1 = (fq_rec_t*) malloc(sizeof(fq_rec_t));
    fq_rec1 = (fq_rec_t*) malloc(sizeof(fq_rec_t));

    barcode_data_t *curr;

    thread_data_t* thread_data = (thread_data_t*)arg;
    int my_line_num;


    /* Get reads, one at a time */

    while(1) {
        // lock reading
        pthread_mutex_lock(thread_data->in_lock);

	//this is equivalent to if(false), which means this block
	//is always skipped, unless when there is an error/end of the file
        if(get_read(fq_rec1, thread_data->params->fq1_fd)) {
	    // sanity check no more reads
            pthread_mutex_unlock(thread_data->in_lock);
	    break;	
	}

        if(thread_data->params->paired > 0) {
            if(get_read(fq_rec2, thread_data->params->fq2_fd)) {
		//error out there becuase if reached the end of the file
		//then we should hit first break, above, since the assumptions
		//that the files of equal length. If issues with R2 only this is an error
                fprintf (stderr, "\n\
                                  \n ERROR: R2 file shorter than R1 file.\
                                  \n Stopping here:\
                                  \n %s\
                                  \n",
                                  fq_rec1->name);
	        //should this be an error?  
	        pthread_mutex_unlock(thread_data->in_lock);
	        exit(1);	
	    }
        }

        // unlock reading
        my_line_num = *(thread_data->line_num);
        *thread_data->line_num += 1;
        pthread_mutex_unlock(thread_data->in_lock);

        int n_crop = 0;

        char *actl_bc = NULL;

        char *fqread1 = NULL;
        char *fqread2 = NULL;

        size_t fq_size = 0;

        fq_size += strlen(fq_rec1->seq);
        fq_size += (strlen(fq_rec1->name)*2);
        fq_size += strlen(fq_rec1->qual);
        fq_size += (strlen(fq_rec1->comment)*2);
        fq_size += 2;// header signs @ and +
        fq_size += 2;//two colons (:)
        fq_size += 4;//cariage returns
        fq_size += 2;//two spaces
        fq_size += 1000;//test

        if(thread_data->params->combine > 0) {
            fq_size += strlen(fq_rec2->seq);
        }

        /* Step 1: Find matching barcode */
        curr = thread_data->curr;
        while(curr) {
            n_crop = chk_bc_mtch(curr->bc, fq_rec1->seq, thread_data->params->mismatch, thread_data->params->max_5prime_crop);
            if(n_crop >= 0) {
                //found matching barcode
                actl_bc = strndup( (fq_rec1->seq)+n_crop, strlen(curr->bc) );
                break;
            }
            curr = curr->next;
        }

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

                const char *actl_umi_idx = (fq_rec1->seq)+strlen(curr->bc)+n_crop;

                if(strlen(actl_umi_idx) < thread_data->params->min_umi_len) {
			//protect by mutex umis_2_short_file
                    fprintf(thread_data->params->umis_2_short_fd, "%s\t%s\t%zu\t%d\n", fq_rec1->name, actl_umi_idx, strlen(actl_umi_idx), thread_data->params->min_umi_len);
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

                get_merged_fqread(&fqread1, fq_rec1, fq_rec2, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);
			//protect by mutex umis_2_short_file
                gzwrite(curr->bcfile1, fqread1, strlen(fqread1));
            }
            else {
                fqread1 = (char*) malloc(fq_size + 1);
                fqread2 = (char*) malloc(fq_size + 1);

                fqread1[0] = '\0';
                fqread2[0] = '\0';

                get_fqread(&fqread1, fq_rec1, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);
                gzwrite(curr->bcfile1, fqread1, strlen(fqread1));

                if(thread_data->params->paired > 0) {
                    get_fqread(&fqread2, fq_rec1, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);
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

            get_fqread(&fqread1, fq_rec1, NULL, NULL, thread_data->params->no_comment, 0);
            gzwrite(thread_data->params->unassigned1_fd, fqread1, strlen(fqread1));
            thread_data->metrics->num_unknown += 1;

            if(thread_data->params->paired > 0) {
                get_fqread(&fqread2, fq_rec2, NULL, NULL, thread_data->params->no_comment, 0);
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
