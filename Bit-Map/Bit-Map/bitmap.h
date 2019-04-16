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
	: m_min(min), m_max(max), m_bit(NULL)	//������ʼ��ʱֱ�Ӹ�ֵ
{
	//����λͼ��ռ�ַ�����С
	int size = (max - min + 8) >> 3;
	m_bit = new unsigned char[size];
	//�����ڴ�ʧ�ܣ������쳣�жϷ���
	if (m_bit == NULL) {
		abort();
	}
	//������λ��0
	memset(m_bit, 0, size);
}

BitMap::~BitMap()
{
	//ʹ��delete[]�ͷ�m_bitȫ���Ĵ洢�ռ䣬deleteֻ���ͷ�m_bit�Ŀռ䣬�����ͷ�ռ�õ���Դ�����ļ���д�ȣ�
	delete[] m_bit;
}

bool BitMap::set(int value)
{
	//�洢�����ֲ���Ԥ�跶Χ��
	if (value < m_min || value > m_max) {
		return(false);
	}
	//���������㿪ʼ���
	value -= m_min;
	//�ҵ�Ҫ�޸ĵ��ֽڣ�һ���ֽڼ�¼8λ������value/8,������value>>3�����ֽ�
	//�ҵ�Ҫ�޸ĵ�bitλ��һ���ֽڼ�¼8λ������value%8---value&0x07,��0x01���ƶ���λ��׼������λ��1
	m_bit[value >> 3] |= (0x01 << (value & 0x07));
	return(true);
}

bool BitMap::clear(int value)
{
	if (value < m_min || value > m_max) {
		return(false);
	}
	value -= m_min;
	//�ҵ��޸ĵ��ֽں;����λ�ã�0x01��λȡ���󽻼�����λ��0
	m_bit[value >> 3] &= ~(0x01 << (value & 0x07));
	return(true);
}

bool BitMap::test(int value)
{
	if (value < m_min || value > m_max) {
		return(false);
	}
	value -= m_min;
	//�ҵ��޸ĵ��ֽں;����λ�ã�0x01��λ�󽻵õ���λ��0��1����λͼ���Ƿ��¼����
	return(m_bit[value >> 3] & (0x01 << (value & 0x07)));
}

#endif


