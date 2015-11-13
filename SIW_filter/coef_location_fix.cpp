//位置修正函数，用于修正相邻区域的影响
//cur: 当前的位置, lo[]偏移的位置
//key偏移的判断添加，len偏移的节点个数
void coef_fix_location_5(int cur, int* lo, int i_0, int i_1, int i_2, int i_3, int i_4, int key){
	if (cur + i_0 >= key || cur + i_0 <= -1)lo[0] = -1 - cur;
	else lo[0] = i_0;
	if (cur + i_1 >= key || cur + i_1 <= -1)lo[1] = -1 - cur;
	else lo[1] = i_1;
	if (cur + i_2 >= key || cur + i_2 <= -1)lo[2] = -1 - cur;
	else lo[2] = i_2;
	if (cur + i_3 >= key || cur + i_3 <= -1)lo[3] = -1 - cur;
	else lo[3] = i_3;
	if (cur + i_4 >= key || cur + i_4 <= -1)lo[4] = -1 - cur;
	else lo[4] = i_4;
}
void coef_fix_location_6(int cur, int* lo, int i_0, int i_1, int i_2, int i_3, int i_4, int i_5, int key){
	if (cur + i_0 >= key || cur + i_0 <= -1)lo[0] = -1 - cur;
	else lo[0] = i_0;
	if (cur + i_1 >= key || cur + i_1 <= -1)lo[1] = -1 - cur;
	else lo[1] = i_1;
	if (cur + i_2 >= key || cur + i_2 <= -1)lo[2] = -1 - cur;
	else lo[2] = i_2;
	if (cur + i_3 >= key || cur + i_3 <= -1)lo[3] = -1 - cur;
	else lo[3] = i_3;
	if (cur + i_4 >= key || cur + i_4 <= -1)lo[4] = -1 - cur;
	else lo[4] = i_4;
	if (cur + i_5 >= key || cur + i_5 <= -1)lo[5] = -1 - cur;
	else lo[5] = i_5;
}
void coef_fix_location_8(int cur, int* lo, int i_0, int i_1, int i_2, int i_3, int i_4, int i_5, int i_6, int i_7, int key){
	if (cur + i_0 >= key || cur + i_0 <= -1)lo[0] = -1 - cur;
	else lo[0] = i_0;
	if (cur + i_1 >= key || cur + i_1 <= -1)lo[1] = -1 - cur;
	else lo[1] = i_1;
	if (cur + i_2 >= key || cur + i_2 <= -1)lo[2] = -1 - cur;
	else lo[2] = i_2;
	if (cur + i_3 >= key || cur + i_3 <= -1)lo[3] = -1 - cur;
	else lo[3] = i_3;
	if (cur + i_4 >= key || cur + i_4 <= -1)lo[4] = -1 - cur;
	else lo[4] = i_4;
	if (cur + i_5 >= key || cur + i_5 <= -1)lo[5] = -1 - cur;
	else lo[5] = i_5;
	if (cur + i_6 >= key || cur + i_6 <= -1)lo[6] = -1 - cur;
	else lo[6] = i_6;
	if (cur + i_7 >= key || cur + i_7 <= -1)lo[7] = -1 - cur;
	else lo[7] = i_7;
}

//switch...case...中由于数组不能初始化失败
//位置修正函数，用于修正相邻区域的影响
//cur: 当前的位置, lo[]偏移的位置，type偏移的类型，
//key偏移的判断添加，len偏移的节点个数
void coef_fix_location(int cur,int* lo, int type, int key, int len){
	if (type == 1){
		for (int i = 0; i != len; i++){
			if (lo[i] + cur > key)
				lo[i] = -1;
		}
	}
	else{
		for (int i = 0; i != len; i++){
			if (lo[i] + cur < key)
				lo[i] = -1;
		}
	}
}