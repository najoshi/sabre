
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

void* demult_runner(void *arg) {

    char fqread1[MAX_READ_SIZE];
    char fqread2[MAX_READ_SIZE];

    fqread1[0] = '\0';
    fqread2[0] = '\0';

    fq_rec_t *fq_rec1 = (fq_rec_t*) malloc(sizeof(fq_rec_t));
    fq_rec_t *fq_rec2 = (fq_rec_t*) malloc(sizeof(fq_rec_t));

    init_fq_rec(fq_rec1);
    init_fq_rec(fq_rec2);

    barcode_data_t *curr;
    curr = NULL;

    thread_data_t* thread_data = (thread_data_t*)arg;
    //int my_line_num;

    /* Get reads, one at a time */

    while(1) {

        // lock reading
        pthread_mutex_lock(thread_data->in_lock);

        //this is equivalent to if(false), which means this block
        //is always skipped, unless when there is an error/end of the file
        if(get_fq_rec(fq_rec1, thread_data->params->fq1_fd)) {
            // sanity check no more reads
            pthread_mutex_unlock(thread_data->in_lock);
            break;
        }

        if(thread_data->params->paired > 0) {
            if(get_fq_rec(fq_rec2, thread_data->params->fq2_fd)) {
                //error out there becuase if reached the end of the file
                //then we should hit first break, above, since the assumptions
                //that the files of equal length. If issues with R2 only this is an error
                fprintf (stderr, "\n\
                        \n ERROR: R2 file shorter than R1 file.\
                        \n Stopping here:\
                        \n %s\
                        \n",
                        fq_rec1->name);
                pthread_mutex_unlock(thread_data->in_lock);
                exit(1);
            }
        }

        // unlock reading
        // TODO this bit of code for ordered fastq files, implement later?
        //my_line_num = *(thread_data->line_num);
        //*thread_data->line_num += 1;
        pthread_mutex_unlock(thread_data->in_lock);

        int n_crop = 0;

        char *actl_bc = NULL;

        /* Step 1: Find matching barcode */
        int got_match = 0;
        curr = thread_data->curr;
        while(curr) {

            for (int i=0; curr->bc[i]; i++) {
                n_crop = chk_bc_mtch(curr->bc[i], fq_rec1->seq, thread_data->params->mismatch, thread_data->params->max_5prime_crop);
                if(n_crop >= 0) {
                    //found matching barcode
                    actl_bc = strndup( (fq_rec1->seq)+n_crop, strlen(curr->bc[i]) );
                    got_match = 1;
                    break;
                }

            }

            if(got_match) {
                break;
            }

            curr = curr->next;
        }

        /* Step 2: Write read out into barcode specific file */
        //TODO this bit of code to keep fastq files ordered as per original fastq files
        //which I don't think that needed? at least not at this stage
        // lock writing
        //while(*(thread_data->out_line_num) != my_line_num) {
        //    pthread_cond_wait(thread_data->cv, thread_data->out_lock);
        //}
        //*thread_data->out_line_num += 1;

        //pthread_cond_broadcast(thread_data->cv);  // Tell everyone it might be their turn!

        char *umi_idx = NULL;

        if(curr != NULL) {
            //for now assume barcode and umi are in R1 read
            if(thread_data->params->umi > 0) {

                const char *actl_umi_idx = (fq_rec1->seq)+strlen(actl_bc)+n_crop;

                if(strlen(actl_umi_idx) < thread_data->params->min_umi_len) {
                    pthread_mutex_lock(thread_data->out_lock);
                    fprintf(thread_data->params->umis_2_short_fd, "%s\t%s\t%zu\t%d\n", fq_rec1->name, actl_umi_idx, strlen(actl_umi_idx), thread_data->params->min_umi_len);
                    pthread_mutex_unlock(thread_data->out_lock);
                    continue;
                }
                else {
                    umi_idx = strdup(actl_umi_idx);
                    umi_idx[thread_data->params->min_umi_len] = '\0';
                }
            }

            if(thread_data->params->combine > 0 && actl_bc != NULL) {
                get_merged_fqread(fqread1, fq_rec1, fq_rec2, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);

                pthread_mutex_lock(thread_data->out_lock);
                fprintf(curr->bcfile1, fqread1);
                pthread_mutex_unlock(thread_data->out_lock);

            }
            else {
                get_fqread(fqread1, fq_rec1, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);

                pthread_mutex_lock(thread_data->out_lock);
                fprintf(curr->bcfile1, fqread1);
                pthread_mutex_unlock(thread_data->out_lock);

                if(thread_data->params->paired > 0) {
                    get_fqread(fqread2, fq_rec1, actl_bc, umi_idx, thread_data->params->no_comment, n_crop);

                    pthread_mutex_lock(thread_data->out_lock);
                    fprintf(curr->bcfile2, fqread2);
                    pthread_mutex_unlock(thread_data->out_lock);

                    //dont need to increment buff_cnt, assuming fq_read1 keeps the right count
                    curr->num_records += 1;
                }
            }
            curr->num_records += 1;
            free(actl_bc);
        }
        else {

            get_fqread(fqread1, fq_rec1, NULL, NULL, thread_data->params->no_comment, 0);

            pthread_mutex_lock(thread_data->out_lock);
            fprintf(thread_data->params->unassigned1_fd, fqread1);
            pthread_mutex_unlock(thread_data->out_lock);

            thread_data->metrics->num_unknown += 1;

            if(thread_data->params->paired > 0) {
                get_fqread(fqread2, fq_rec2, NULL, NULL, thread_data->params->no_comment, 0);

                pthread_mutex_lock(thread_data->out_lock);
                fprintf(thread_data->params->unassigned2_fd, fqread2);
                pthread_mutex_unlock(thread_data->out_lock);

                thread_data->metrics->num_unknown += 1;
            }

        }

        thread_data->metrics->total += 2;
    }

    free(fq_rec1);
    free(fq_rec2);
    //free(thread_data); according to valgrind report this line isn't needed since no errors given out..
    return NULL;
}
