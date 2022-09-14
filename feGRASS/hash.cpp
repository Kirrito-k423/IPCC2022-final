#include "global.h"

unsigned int _hash(unsigned int x){
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

/**
 * 将元素x放入哈希表中，冲突时尝试放入下一个位置。
 * 成功放入返回1，失败（哈希表满了）返回0
 */
int hash_put(simple_map *hash_table_, int start_offset,unsigned int key, int value){
	simple_map *hash_table = hash_table_ + start_offset * HASH_SIZE;
    unsigned int idx = _hash(key) % HASH_SIZE;
    if(hash_table[idx].key==EMPTY){
        hash_table[idx].key = key;
		hash_table[idx].value = value;
        return 1;
    }
    int cnt=0;
    do{
		DEBUG_PRINT("hash_put_hit ");
        idx = (idx + 1) % HASH_SIZE;
        if(hash_table[idx].key==EMPTY){
			hash_table[idx].key = key;
			hash_table[idx].value = value;
			DEBUG_PRINT("\n");
            return 1;
        }
        cnt++;
    }while(cnt<HASH_SIZE);
    return 0;
}


/**
 * 判断元素x是否在哈希表中
 * 在返回1，不在返回0
 */
int hash_get_value(simple_map *hash_table_, int start_offset,unsigned int key){
	simple_map *hash_table = hash_table_ + start_offset * HASH_SIZE;
    unsigned int idx = _hash(key) % HASH_SIZE;
    simple_map tmp = hash_table[idx];
    if(tmp.key==EMPTY){
        return -1;
    }else if(tmp.key==key){
        return tmp.value;
    }else{
        int cnt=0;
        do{
			// DEBUG_PRINT("hash_get_hit ");
            idx = (idx + 1) % HASH_SIZE;
            if(hash_table[idx].key==key){
                return hash_table[idx].value;
            }
            cnt++;
        }while(cnt<HASH_SIZE);
        return -1;
    }
}