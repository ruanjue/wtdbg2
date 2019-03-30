#ifndef __TEXT_PLOT_RJ_H
#define __TEXT_PLOT_RJ_H

#include "list.h"

static inline char* barplot_txt_u4_simple(u4i w, u4i h, u4i *vals, u4i size, u4i max_val){
	char *g;
	double hdiv;
	u4i wdiv;
	u4i i, j, x, y;
	g = malloc((w + 1) * h + 1);
	memset(g, ' ', (w + 1) * h);
	for(i=1;i<=h;i++){
		g[(w + 1) * i - 1] = '\n';
	}
	g[(w + 1) * h] = 0;
	if(max_val == 0){
		max_val = 1;
		for(i=0;i<size;i++){
			if(vals[i] > max_val) max_val = vals[i];
		}
	}
	hdiv = 1.0 * max_val / h;
	wdiv = size / w;
	if(wdiv == 0) wdiv = 1;
	for(i=x=0;i<size;i+=wdiv,x++){
		y = UInt(vals[i] / hdiv);
		if(y == 0 && vals[i]) y = 1;
		if(y > h) y = h;
		for(j=1;j<=y;j++){
			g[(h - j) * (w + 1) + x] = '|';
		}
	}
	return g;
}

static inline char* barplot_txt_u8_simple(u4i w, u4i h, u8i *vals, u4i size, u8i max_val){
	char *g;
	double hdiv;
	u8i wdiv;
	u4i i, j, x, y;
	g = malloc((w + 1) * h + 1);
	memset(g, ' ', (w + 1) * h);
	for(i=1;i<=h;i++){
		g[(w + 1) * i - 1] = '\n';
	}
	g[(w + 1) * h] = 0;
	if(max_val == 0){
		max_val = 1;
		for(i=0;i<size;i++){
			if(vals[i] > max_val) max_val = vals[i];
		}
	}
	hdiv = 1.0 * max_val / h;
	wdiv = size / w;
	if(wdiv == 0) wdiv = 1;
	for(i=x=0;i<size;i+=wdiv,x++){
		y = UInt(vals[i] / hdiv);
		if(y == 0 && vals[i]) y = 1;
		if(y > h) y = h;
		for(j=1;j<=y;j++){
			g[(h - j) * (w + 1) + x] = '|';
		}
	}
	return g;
}

#endif
