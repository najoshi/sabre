#include "utils.h"

// https://stackoverflow.com/questions/2336242/recursive-mkdir-system-call-on-unix/11425692
// https://stackoverflow.com/questions/7430248/creating-a-new-directory-in-c
const char * _mkdir(const char *file_path) {

    if(!file_path) {
        fprintf (stderr,
                "ERROR: This shouldn't happend, file path == %s\n",
                file_path);
        exit(EXIT_FAILURE);
    }

    // return straigth away if a file_path is not nested file path
    if(strstr(file_path, "/") == NULL) {
        return file_path;
    }
    //TODO check if directory already exists or not
    //struct stat st = {0};
    //
    //if (stat(dir, &st) == -1) {
    //    mkdir(tmp, S_IRWXU);
    //}

    //char tmp[PATH_MAX]; // can't get this to work..
    char tmp[256];
    char *p = NULL;
    char *dirc = strdup(file_path);
    size_t len;

    const char *dir = dirname(dirc);
    snprintf(tmp, sizeof(tmp),"%s", dir);
    len = strlen(tmp);

    if(tmp[len - 1] == '/') {
        tmp[len - 1] = 0;
    }

    for(p = tmp + 1; *p; p++) {
        if(*p == '/') {
            *p = 0;
            mkdir(tmp, S_IRWXU);
            *p = '/';
        }
    }

    mkdir(tmp, S_IRWXU);
    free(dirc);

    return file_path;
}

//NOTE retuns zero on success
//strcmp can be used for sorting, returns pos, zero, neg
//BUT this new implementation can't be used as such just FYI
match_ret_t chk_bc_mtch(const char *orig_bc, const char *orig_read, size_t max_mismatch, int max_5prime_crop) {
    int orig_read_len = strlen(orig_read);
    int orig_bc_len = strlen(orig_bc);
    int n_crop = 0;
    match_ret_t ret = { -1,-1 };
    if(orig_bc_len > orig_read_len) {
        fprintf (stderr,
                "WARNING: Length of the barcode %d is greater than length of the reads %d.",
                orig_bc_len, orig_read_len);
        return ret;
    }

    while(n_crop <= max_5prime_crop) {

        if(n_crop > orig_read_len) {
            return ret;
        }

        int mismatch = 0;
        char u1, u2;
        const char *bc = orig_bc;
        const char *read = orig_read+n_crop;
        int bc_len = orig_bc_len;

        while (bc_len-- > 0) {
            u1 = *bc++;
            u2 = *read++;

            if (u1 != u2) {
                mismatch++;
                if (mismatch > max_mismatch) {
                    break;
                }
            }

            if (u1 == '\0' || u2 == '\0') {
                ret.cropped = n_crop;
                ret.mismatches = mismatch;
                return ret;
            }
        }

        if(mismatch <= max_mismatch) {
            ret.cropped = n_crop;
            ret.mismatches = mismatch;
            return ret;
        }

        n_crop++;
    }
    //this is in the case of error
    return ret;
}

// https://stackoverflow.com/questions/21880730/c-what-is-the-best-and-fastest-way-to-concatenate-strings
//TODO this is a fastq mystrcat function, that returns a pointer to the end of the string
void get_fqread(char *fqread, fq_rec_t *fq_rec, char *barcode, char *umi_idx, int no_comment, int n_crop) {

    fqread[0] = '\0';

    if(n_crop < 0) {
        fprintf(stderr,
                "ERROR: n_crop set to %d. This can't happend\n",
                n_crop);
        exit(EXIT_FAILURE);
    }

    //@READNAME:BACRCODE:UMI
    //1st line
    strcat(fqread, fq_rec->name);
    //TODO later can have conditional here depending on the the structure and/or BARCODE/UMI
    if(barcode) {
        strcat(fqread, ":");
        strcat(fqread, barcode);
    }
    else if(!barcode) {
        barcode = "";
    }

    if(umi_idx) {
        strcat(fqread, ":");
        strcat(fqread, umi_idx);
    }

    if(fq_rec->comment && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fq_rec->comment);
    }
    strcat(fqread, "\n");

    //2nd line
    strcat(fqread, (fq_rec->seq)+strlen(barcode)+n_crop);
    strcat(fqread, "\n");

    //3rd line
    strcat(fqread, "+");
    strcat(fqread, fq_rec->name);
    if(fq_rec->comment && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fq_rec->comment);
    }
    strcat(fqread, "\n");

    //4th line
    strcat(fqread, (fq_rec->qual)+strlen(barcode)+n_crop);
    strcat(fqread, "\n");
}

void get_merged_fqread(char *fqread, fq_rec_t *fq_rec1, fq_rec_t *fq_rec2, char *barcode, char *umi_idx, int no_comment, int n_crop) {
    fqread[0] = '\0';
    //@READNAME:BACRCODE:UMI
    //1st line
    strcat(fqread, fq_rec1->name);
    //TODO later can have conditional here depending on the the structure and/or BARCODE/UMI
    if(barcode) {
        strcat(fqread, ":");
        strcat(fqread, barcode);
    }

    if(umi_idx) {
        strcat(fqread, ":");
        strcat(fqread, umi_idx);
    }

    if(fq_rec1->comment && no_comment == -1) {
        strcat(fqread, " ");
        strcat(fqread, fq_rec1->comment);
    }
    strcat(fqread, "\n");

    //2nd line
    strcat(fqread, fq_rec2->seq);
    strcat(fqread, "\n");

    //3rd line
    strcat(fqread, "+");
    strcat(fqread, "");
    //strcat(fqread, fq_rec2->name);
    //if(fq_rec2->comment && no_comment == -1) {
    //    strcat(fqread, " ");
    //    strcat(fqread, fq_rec2->comment);
    //}
    strcat(fqread, "\n");

    //4th line
    strcat(fqread, (fq_rec2->qual));
    strcat(fqread, "\n");
}

void get_bc_fn(char **bcout_fn, char *s_name, char *barcode, int read_type) {

    if(strlen(s_name) > MAX_FILENAME_LENGTH) {
        fprintf (stderr,
                "ERROR: Too many characters in your sample name; %s:%zd \n",
                s_name, strlen(s_name));
        exit(EXIT_FAILURE);
    }

    strcat(*bcout_fn, s_name);
    strcat(*bcout_fn, "_");
    strcat(*bcout_fn, barcode);

    if(read_type == 1) {
        strcat(*bcout_fn, "_R1.fastq");
    }
    else if(read_type == 2) {
        strcat(*bcout_fn, "_R2.fastq");
    }
    else {
        fprintf (stderr,
                "ERROR: This shouldn't happen, wrong read type was passed through -> %d\n",
                read_type);
        exit(EXIT_FAILURE);
    }
}

void set_default_params(param_t *params) {
    params->mismatch = 0;
    params->combine = -1;
    params->umi = -1;
    params->paired = -1;
    params->min_umi_len = 0;
    params->max_5prime_crop = 0;
    params->no_comment = -1;
}

void params_destroy(param_t *params) {
    gzclose(params->fq1_fd);
    gzclose(params->fq2_fd);
    fclose(params->unassigned1_fd);
    fclose(params->unassigned2_fd);
    fclose(params->umis_2_short_fd);
}
