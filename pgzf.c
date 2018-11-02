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
	"Version: 1.0\n"
	"Usage: pgzf [options] file1 [file2 ...]\n"
	"Options:\n"
	" -d          Decompress mode\n"
	" -t <int>    Number of threads, [8]\n"
	" -f          Force to overwrite\n"
	" -o <string> Output file name\n"
	" -x          Delete input files after done, need to specify `-o` or/and `-f`\n"
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
	char *outf;
	FILE *in, *out;
	void *buff;
	u4i bufsize, nbyte;
	int c, rw, ncpu, level, overwrite, del;
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
			case 'V': fprintf(stdout, "pgzf 1.0\n"); return 0;
			default: return usage(1);
		}
	}
	if(del){
		if(outf == NULL && overwrite == 0){
			if(optind < argc){
				fprintf(stderr, " ** WARNNING: won't delete input files. To force delete input files, please specify -o or/and -f\n");
			}
			del = 0;
		}
	}
	if(outf){
		for(c=optind;c<argc;c++){
			if(strcmp(outf, argv[c]) == 0){
				fprintf(stderr, " ** ERROR: The same file in INPUT and OUTPUT, '%s'\n", outf);
				return 1;
			}
		}
		out = open_file_for_write(outf, NULL, overwrite);
	} else {
		out = stdout;
	}
	buff = malloc(bufsize);
	if(rw == PGZF_MODE_R){
		do {
			if(optind < argc){
				in = open_file_for_read(argv[optind], NULL);
			} else {
				in = stdin;
			}
			pz = open_pgzf_reader(in, bufsize, ncpu);
			while((nbyte = read_pgzf(pz, buff, bufsize))){
				fwrite(buff, 1, nbyte, out);
			}
			close_pgzf(pz);
			if(in != stdin){
				fclose(in);
				if(del){
					unlink(argv[optind]);
				}
			}
			optind ++;
		} while(optind < argc);
	} else {
		pz = open_pgzf_writer(out, bufsize, ncpu, level);
		do {
			if(optind < argc){
				in = open_file_for_read(argv[optind], NULL);
			} else {
				in = stdin;
			}
			while((nbyte = fread(buff, 1, bufsize, in))){
				write_pgzf(pz, buff, nbyte);
			}
			if(in != stdin){
				fclose(in);
				if(del){
					unlink(argv[optind]);
				}
			}
			optind ++;
		} while(optind < argc);
		close_pgzf(pz);
	}
	if(out != stdout) fclose(out);
	free(buff);
	return 0;
}
