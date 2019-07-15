/*
 *
 * Copyright (c) 2018, Jue Ruan <ruanjue@gmail.com>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "pgzf.h"

int usage(int ret){
	fprintf(stdout,
	"PGZF: Parallel gzip file IO\n"
	"Author: Jue Ruan <ruanjue@gmail.com>\n"
	"Version: 1.1\n"
	"Usage: pgzf [options] file1 [file2 ...]\n"
	"Options:\n"
	" -d          Decompress mode\n"
	" -t <int>    Number of threads, [8]\n"
	" -f          Force to overwrite\n"
	" -o <string> Output file name, support directory\n"
	" -x          Delete input files after done\n"
	" -b <int>    Block size in MB, 1 ~ 256 [16]\n"
	" -l <int>    Compress level, 1-9, see gzip, [6]\n"
	" -h          Show this document\n"
	" -V          Print version information and exit\n"
	"\n"
	"File format:\n"
	" PGZF fellows standard GZIP format (rfc1952), and is blocked compressed.\n"
	" It defines two TAGs in each GZIP header, ZS: block size, ZX: random access index.\n"
	" Program pgzf can decompress .pgzf and .gz files. When decompressing .gz files,\n"
	" pgzf is in fact a buffered gzip reader. Also, .pgzf files can be decompressed\n"
	" by program gzip.\n"
	"\n"
	"In plan to support random access\n"
	);
	return ret;
}

int main(int argc, char **argv){
	PGZF *pz;
	char *outf, *ftag;
	FILE *in, *out;
	void *buff;
	u4i bufsize, nbyte;
	int c, rw, ncpu, level, overwrite, del, is_dir;
	rw = PGZF_MODE_W;
	ncpu = 8;
	bufsize = PGZF_DEFAULT_BUFF_SIZE;
	level = 6;
	overwrite = 0;
	del = 0;
	outf = NULL;
	while((c = getopt(argc, argv, "hdxft:b:l:o:V")) != -1){
		switch(c){
			case 'h': return usage(0);
			case 'd': rw = PGZF_MODE_R; break;
			case 't': ncpu = atoi(optarg); break;
			case 'b': bufsize = (atol(optarg) << 20); break;
			case 'l': level = atoi(optarg); break;
			case 'f': overwrite = 1; break;
			case 'o': outf = optarg; break;
			case 'x': del = 1; break;
			case 'V': fprintf(stdout, "pgzf 1.1\n"); return 0;
			default: return usage(1);
		}
	}
	if(optind == argc){
		return usage(1);
	}
	if(0 && del){
		if(outf == NULL && overwrite == 0){
			if(optind < argc){
				fprintf(stderr, " ** WARNNING: won't delete input files. To force delete input files, please specify -o or/and -f\n");
			}
			del = 0;
		}
	}
	is_dir = 0;
	out = NULL;
	if(outf){
		if(file_exists(outf)){
			if(overwrite == 0){
				fprintf(stderr, " ** ERROR: '%s' exists\n", outf);
				return 1;
			} else {
				for(c=optind;c<argc;c++){
					if(strcmp(outf, argv[c]) == 0){
						fprintf(stderr, " ** ERROR: The same file in INPUT and OUTPUT, '%s'\n", outf);
						return 1;
					}
				}
			}
			out = open_file_for_write(outf, NULL, overwrite);
		} else if(dir_exists(outf)){
			is_dir = 1;
		} else {
			out = open_file_for_write(outf, NULL, overwrite);
		}
	}
	buff = malloc(bufsize);
	if(rw == PGZF_MODE_R){
		if(outf == NULL || is_dir){
			for(c=optind;c<argc;c++){
				if(strlen(argv[c]) < 4 || strcasecmp(argv[c] + strlen(argv[c]) - 3, ".gz")){
					fprintf(stderr, " ** ERROR: cannot auto generate output file name for '%s'\n", argv[c]);
					return 1;
				} else if(is_dir){
					char *rtag;
					rtag = relative_filename(argv[c]);
					rtag[strlen(rtag) - 3] = 0;
					ftag = malloc(strlen(outf) + 1 + strlen(rtag) + 1);
					sprintf(ftag, "%s/%s", outf, rtag);
					free(rtag);
					if(overwrite == 0 && file_exists(ftag)){
						fprintf(stderr, " ** ERROR: '%s' exists\n", ftag);
						return 1;
					}
					free(ftag);
				} else {
					ftag = strdup(argv[optind]);
					ftag[strlen(ftag) - 3] = 0;
					if(overwrite == 0 && file_exists(ftag)){
						fprintf(stderr, " ** ERROR: '%s' exists\n", ftag);
						return 1;
					}
					free(ftag);
				}
			}
		}
		do {
			in = open_file_for_read(argv[optind], NULL);
			if(outf == NULL){
				ftag = strdup(argv[optind]);
				ftag[strlen(ftag) - 3] = 0;
				out = open_file_for_write(ftag, NULL, overwrite);
				free(ftag);
			} else if(is_dir){
				char *rtag;
				rtag = relative_filename(argv[optind]);
				rtag[strlen(rtag) - 3] = 0;
				ftag = malloc(strlen(outf) + 1 + strlen(rtag) + 1);
				sprintf(ftag, "%s/%s", outf, rtag);
				free(rtag);
				out = open_file_for_write(ftag, NULL, overwrite);
				free(ftag);
			}
			pz = open_pgzf_reader(in, bufsize, ncpu);
			while((nbyte = read_pgzf(pz, buff, bufsize))){
				fwrite(buff, 1, nbyte, out);
			}
			if(pz->error){
				fprintf(stderr, " ** ERROR: error code (%d)'\n", pz->error);
				return 1;
			}
			close_pgzf(pz);
			if(in != stdin){
				fclose(in);
				if(del){
					unlink(argv[optind]);
				}
			}
			optind ++;
			if(outf == NULL || is_dir){
				fclose(out);
			}
		} while(optind < argc);
	} else {
		if(outf && !is_dir){
			pz = open_pgzf_writer(out, bufsize, ncpu, level);
		} else {
			pz = NULL;
			for(c=optind;c<argc;c++){
				if(strlen(argv[c]) >= 4 && strcasecmp(argv[c] + strlen(argv[c]) - 3, ".gz") == 0){
					fprintf(stderr, " ** ERROR: file seems already compressed '%s'\n", argv[c]);
					return 1;
				} else if(strcmp(argv[c], "-") == 0){
					fprintf(stderr, " ** ERROR: Please specify output file when read from STDIN '%s'\n", argv[c]);
					return 1;
				} else if(is_dir){
					char *rtag;
					rtag = relative_filename(argv[c]);
					ftag = malloc(strlen(outf) + 1 + strlen(rtag) + 3 + 1);
					sprintf(ftag, "%s/%s.gz", outf, rtag);
					free(rtag);
					if(overwrite == 0 && file_exists(ftag)){
						fprintf(stderr, " ** ERROR: '%s' exists\n", ftag);
						return 1;
					}
					free(ftag);
				} else {
					ftag = malloc(strlen(argv[c]) + 4);
					sprintf(ftag, "%s.gz", argv[c]);
					if(overwrite == 0 && file_exists(ftag)){
						fprintf(stderr, " ** ERROR: '%s' exists\n", ftag);
						return 1;
					}
					free(ftag);
				}
			}
		}
		do {
			if(outf == NULL){
				ftag = malloc(strlen(argv[optind]) + 4);
				sprintf(ftag, "%s.gz", argv[optind]);
				out = open_file_for_write(ftag, NULL, overwrite);
				pz = open_pgzf_writer(out, bufsize, ncpu, level);
				free(ftag);
			} else if(is_dir){
				char *rtag;
				rtag = relative_filename(argv[optind]);
				ftag = malloc(strlen(outf) + 1 + strlen(rtag) + 3 + 1);
				sprintf(ftag, "%s/%s.gz", outf, rtag);
				free(rtag);
				out = open_file_for_write(ftag, NULL, overwrite);
				pz = open_pgzf_writer(out, bufsize, ncpu, level);
				free(ftag);
			}
			in = open_file_for_read(argv[optind], NULL);
			while((nbyte = fread(buff, 1, bufsize, in))){
				write_pgzf(pz, buff, nbyte);
			}
			if(in != stdin){
				fclose(in);
				if(del){
					unlink(argv[optind]);
				}
			}
			if(outf == NULL || is_dir){
				close_pgzf(pz);
				fclose(out);
			}
			optind ++;
		} while(optind < argc);
		if(outf && !is_dir){
			close_pgzf(pz);
		}
	}
	if(outf && !is_dir) fclose(out);
	free(buff);
	return 0;
}
