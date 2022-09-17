#include "global.h"

static int *bitmap;
static int m;   //m bits
static int n;   //n element
static int bit_patterns[32];
// std::hash<int64_t> hasher;

unsigned int my_hash(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

static inline int hash_(int u, int v){
    // int h = hasher(((int64_t)u<<32)|v);
    // int h = my_hash(u * v);
    int h = u^v;    //key为(u,v)，简单的hash
    return h;
}

/**
 * bloom filter
 * n为预计插入元素个数
 * 分配n*32位bit，之后插入n个元素，使用1个hash函数，假阳性率约为0.03
*/
void bloom_init(int n_){
    n = n_;
    m = n*32;
    bitmap = (int *)malloc(n*sizeof(int));
    memset(bitmap, 0, n*sizeof(int));

    for(int i=0; i<32; i++){
        bit_patterns[i] = 1<<i;
    }
}

void bloom_insert(int u, int v){
    int h = hash_(u, v);
    bitmap[(h>>5)%n] |= bit_patterns[(h&0x1f)];  //插入元素
}

bool bloom_contain(int u, int v){
    int h = hash_(u, v);
    return bitmap[(h>>5)%n] & bit_patterns[(h&0x1f)];
}

void bloom_free(){
    free(bitmap);
}

void test_bitmap() {
    int hit, total;
    total = m;
    hit = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 32; j++) {
            if (bitmap[i] & bit_patterns[j]) {
                hit++;
            }
        }
    }
    printf("one bit: %d, total: %d, rate: %f\n", hit, total, (double)hit / total);
}