/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-27 15:58:17
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-08-28 00:18:43
 */
#ifndef _GLOBAL_H
#define _GLOBAL_H
// #include <omp.h>

#ifdef TIME
#define TIME_PRINT(fmt, args...)    fprintf(stderr, fmt, ## args)
#else
#define TIME_PRINT(fmt, args...)    /* Don't do anything in release builds */
#endif

#ifdef DEBUG
#define PIntWithName(var)  printf("%s\t%d\n", #var, var) 
#define DEBUG_PRINT(fmt, args...)    fprintf(stderr, fmt, ## args)
#define OMP_TIME_PRINT(fmt, args...)    fprintf(stderr, fmt, ## args)
#else
#define PIntWithName(var)    /* Don't do anything in release builds */
#define DEBUG_PRINT(fmt, args...)  
#define OMP_TIME_PRINT(fmt, args...)
#endif

#endif