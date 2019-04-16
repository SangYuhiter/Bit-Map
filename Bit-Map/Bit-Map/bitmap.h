#pragma once

#ifndef BITMAP_H
#define BITMAP_H

#include <cstdlib>
#include <cstring>

class BitMap
{
public:
	BitMap(int min, int max);
	~BitMap();
	bool set(int value);
	bool clear(int value);
	bool test(int value);
private:
	BitMap(const BitMap &);
	BitMap & operator = (const BitMap &);
private:
	int             m_min;
	int             m_max;
	unsigned char * m_bit;
};

BitMap::BitMap(int min, int max)
	: m_min(min), m_max(max), m_bit(NULL)	//函数初始化时直接赋值
{
	//计算位图所占字符数大小
	int size = (max - min + 8) >> 3;
	m_bit = new unsigned char[size];
	//申请内存失败，调用异常中断返回
	if (m_bit == NULL) {
		abort();
	}
	//将比特位置0
	memset(m_bit, 0, size);
}

BitMap::~BitMap()
{
	//使用delete[]释放m_bit全部的存储空间，delete只能释放m_bit的空间，不能释放占用的资源（如文件读写等）
	delete[] m_bit;
}

bool BitMap::set(int value)
{
	//存储的数字不在预设范围内
	if (value < m_min || value > m_max) {
		return(false);
	}
	//做减法从零开始标记
	value -= m_min;
	//找到要修改的字节，一个字节记录8位数，即value/8,所以是value>>3处的字节
	//找到要修改的bit位，一个字节记录8位数，即value%8---value&0x07,将0x01左移多少位即准备将该位置1
	m_bit[value >> 3] |= (0x01 << (value & 0x07));
	return(true);
}

bool BitMap::clear(int value)
{
	if (value < m_min || value > m_max) {
		return(false);
	}
	value -= m_min;
	//找到修改的字节和具体的位置，0x01移位取反求交即将该位置0
	m_bit[value >> 3] &= ~(0x01 << (value & 0x07));
	return(true);
}

bool BitMap::test(int value)
{
	if (value < m_min || value > m_max) {
		return(false);
	}
	value -= m_min;
	//找到修改的字节和具体的位置，0x01移位求交得到该位，0，1即是位图中是否记录该数
	return(m_bit[value >> 3] & (0x01 << (value & 0x07)));
}

#endif


