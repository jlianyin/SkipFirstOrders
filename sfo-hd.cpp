/* 
 * sfo-hd.cpp - Computes Hilbert space-filling curve coordinates from integer index, and vice versa.  
 * 
 * Author:      Lianyin Jia
 *              Faculty of Information Engineering and Automation
 *              Kunming University of Science & Technology               
 * Date:        October 20 2020
 * Copyright (c) 2020-2022, Kunming University of Science & Technology
 */

#include <stdio.h> 
typedef unsigned long long bitmask_t;
typedef unsigned long halfmask_t;

#define MAX(x,y) (x>y?x:y)
//get the value of a certain bit by a position i
#define getOneBitByPos(X,i) ((X>>i) & 0x01) 
//get 2 continuous bits by a position starting from i
#define getTwoBitByPos(X,i) (X>>(i-1)& 3UL) 

//generate a number with k 1s
#define ones(T,k) ((((T)2) << (k-1)) - 1) 

//state view for encoding
char arKey[4][2][2] = { 0,1,3,2,0,3,1,2,2,3,1,0,2,1,3,0 };
char arType[4][2][2] = { 1,0,3,0,0,2,1,1,2,1,2,3,3,3,0,2 };
//state view for decoding
char invKey[4][4] ={0,1,3,2,0,2,3,1,3,2,0,1,3,1,0,2};
char invType[4][4] ={1,0,0,3,0,1,1,2,3,2,2,1,2,3,3,0};

//Look up table for first m orders
bitmask_t PartialKeys[2][2][32]={0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,
	0x1,0x5,0x15,0x55,0x155,0x555,0x1555,0x5555,0x15555,0x55555,0x155555,0x555555,0x1555555,0x5555555,0x15555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,0x55555555,
	0x3,0xf,0x3f,0xff,0x3ff,0xfff,0x3fff,0xffff,0x3ffff,0xfffff,0x3fffff,0xffffff,0x3ffffff,0xfffffff,0x3fffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,0xffffffff,
	0x2,0xa,0x2a,0xaa,0x2aa,0xaaa,0x2aaa,0xaaaa,0x2aaaa,0xaaaaa,0x2aaaaa,0xaaaaaa,0x2aaaaaa,0xaaaaaaa,0x2aaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa,0xaaaaaaaa};
char PartialStatus[2][2][2]={0,1,0,0,0,3,0,0};
 
//Find First Bit for a 32-bit number
int msb32_idx(halfmask_t n)
{
	int b = 0;
	if (!n) return -1;  
#define step(x) if (n >= ((halfmask_t)1) << x) b += x, n >>= x
	step(16); step(8); step(4); step(2); step(1);
#undef step
	return b;
}
//Find First Bit for a 64-bit number
int msb64_idx(bitmask_t n)
{
	int b = 0;
	if (!n) return -1;  
#define step(x) if (n >= ((bitmask_t)1) << x) b += x, n >>= x
	step(32); step(16); step(8); step(4); step(2); step(1);
#undef step
	return b;
}

/*****************************************************************
 * sfo_he
 * 
 * Convert coordinates of a point on a Hilbert curve to its index.
 * Inputs:
 *  GridX:      X coordinates.
 *  GridY:      Y coordinate.
 *  k:      the number of orders.
 * Outputs:
 *  index:      Output index value.
 */
bitmask_t sfo_he(halfmask_t GridX, halfmask_t GridY, int k)
{
	unsigned bitX = 0, bitY = 0;	
	halfmask_t tempX=GridX, tempY=GridY;
	bitmask_t resKey = 0, A;
	int startPos, m;
	unsigned nType=0;
	halfmask_t mask = ones(halfmask_t,k),mask1 = 1<<k-1;
	
	if(GridX>=mask1){ bitX = 1; tempX=~tempX & mask; }
	if(GridY>=mask1){ bitY = 1; tempY=~tempY & mask; }
	startPos = msb32_idx(MAX(tempX,tempY));
	m=k-startPos-1;
	A = ((1<<2*m)-1)/3;
	if(bitX==1) 
	{
		if(bitY==1) 
			resKey = 2 * A;
		else
		{			
			resKey = 3*A; 
			nType = m%2==0?0:3;
		}		
	}
	else 
	{
		if(bitY==1)		
			resKey = A; 		
		else
			nType = m%2==0?0:1;
	}
	for (int i = startPos; i >= 0; i--)
	{ 
		bitX = getOneBitByPos(GridX,i);
		bitY = getOneBitByPos(GridY,i);
		resKey = (resKey << 2) | arKey[nType][bitX][bitY]; 
		nType = arType[nType][bitX][bitY];
	}	 
	return resKey; 
}

/*****************************************************************
 * sfo_lut_he
 * 
 * Convert coordinates of a point on a Hilbert curve to its index with an additional look up table.
 * Inputs:
 *  GridX:      X coordinates.
 *  GridY:      Y coordinate.
 *  k:          The number of orders.
 * Outputs:
 *  index:      Output index value.
 */
bitmask_t sfo_lut_he(halfmask_t GridX, halfmask_t GridY, int k)
{
	unsigned bitX = 0, bitY = 0;
	halfmask_t tempX=GridX,tempY=GridY;
	halfmask_t mask = (1<<k)-1,mask1 = 1<<k-1; 
	int startPos,m;
	bitmask_t resKey;

	if(GridX>=mask1){ bitX = 1; tempX=~tempX & mask; }
	if(GridY>=mask1){ bitY = 1; tempY=~tempY & mask; }
	startPos = msb32_idx(MAX(tempX,tempY));
	m=k-startPos-1;
	resKey = PartialKeys[bitX][bitY][m-1]; 
	unsigned nType = PartialStatus[bitX][bitY][m%2]; 
	for (int i = startPos; i >= 0; i--)
	{  
		bitX = getOneBitByPos(GridX,i);
		bitY = getOneBitByPos(GridY,i);
		resKey = (resKey << 2) | arKey[nType][bitX][bitY];
		nType = arType[nType][bitX][bitY];
	}	 
	return resKey;
}

/*****************************************************************
 * sfo_de
 * 
 * Convert an index into a Hilbert curve to a set of coordinates.
 * Inputs:
 *  index:      The index.
 *  GridX:      X coordinates.
 *  GridY:      Y coordinates.
 *  k:          The number of orders.
 * Outputs:
 *  GridX and GridY:      The computed coordinates.
 */
void sfo_de(bitmask_t index,halfmask_t &GridX, halfmask_t &GridY, int k)
{ 
	unsigned bitIndx = getTwoBitByPos(index,2*k-1);
	bitmask_t resKey = 0; 
	int startPos,order; 
	unsigned nType=0;
	bitmask_t mask = ones(bitmask_t,2*k);
	char posKey;
	
	switch(bitIndx)
	{
	case 0:		
		startPos = msb64_idx(index);
		order = k - startPos/2 -1;
		nType = order%2;
		break;
	case 1:
		startPos = msb64_idx(~(index ^ index<<1)&mask);
		order = k - (startPos-1)/2 -1; 
		GridY = ones(halfmask_t,order); 
		break;
	case 2:
		startPos = msb64_idx(~(index ^ index>>1)&mask);
		order = (2*k - startPos-1)/2; 
		GridX = GridY = ones(halfmask_t,order); 
		break;
	case 3:
		startPos = msb64_idx(~index & mask);
		order = (2*k - startPos-1)/2; 
		GridX = ones(halfmask_t,order); 		
		nType = order%2==0?0:3;
		break;		 
	}

	for (int i = k-order; i>0; i--)
	{		
		bitIndx = getTwoBitByPos(index,2*i-1);  		
		posKey = invKey[nType][bitIndx]; 
		GridY = GridY <<1 | posKey & 0x1;
		GridX = GridX <<1 | posKey>>1 & 0x1;		
		nType = invType[nType][bitIndx];
	} 
}

int main()
{
	bitmask_t key = sfo_he(3,7,4);
	bitmask_t key1 =sfo_lut_he(3,7,4);

	halfmask_t gridX=0, gridY=0;
	printf("%lld,%lld\n",key,key1);
	 
	sfo_de(key,gridX,gridY,4);
	printf("%d,%d\n",gridX,gridY); 
	
	return 0;

}
