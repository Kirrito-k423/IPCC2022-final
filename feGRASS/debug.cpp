/*
 * @Descripttion: 
 * @version: 
 * @Author: Shaojie Tan
 * @Date: 2022-08-31 22:11:56
 * @LastEditors: Shaojie Tan
 * @LastEditTime: 2022-08-31 22:37:37
 */
#include "global.h"

void printStack(string name, stack<int> toPrint){
    DEBUG_PRINT("%s:\n",name.c_str());
    int i=0;
    while(!toPrint.empty()){
        i++;
        if(i%13==12){
            DEBUG_PRINT("%d\n",toPrint.top());
        }else
            DEBUG_PRINT("%d\t",toPrint.top());
        toPrint.pop();
    }
    DEBUG_PRINT("\n");
}

void print_M1_Array(string name,int * toPrint){
    DEBUG_PRINT("%s:\n",name.c_str());
    for(int i = 0; i < M+1; i++){
        if(i%13==12){
            DEBUG_PRINT("%d\n",toPrint[i]);
        }else{
            DEBUG_PRINT("%d\t",toPrint[i]);
        }
    }
    DEBUG_PRINT("\n");

}