#pragma once
/*
 *********************************************************************
 *                                                                   *
 *                           Open Bloom Filter                       *
 *                                                                   *
 * Author: Arash Partow - 2000                                       *
 * URL: http://www.partow.net                                        *
 * URL: http://www.partow.net/programming/hashfunctions/index.html   *
 *                                                                   *
 * Copyright notice:                                                 *
 * Free use of the Open Bloom Filter Library is permitted under the  *
 * guidelines and in accordance with the MIT License.                *
 * http://www.opensource.org/licenses/MIT                            *
 *                                                                   *
 *********************************************************************
*/


#ifndef INCLUDE_BLOOM_FILTER_HPP
#define INCLUDE_BLOOM_FILTER_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <string>
#include <vector>


static const std::size_t bits_per_char = 0x08;    // 8 bits in 1 char(unsigned)

static const unsigned char bit_mask[bits_per_char] = {
													   0x01,  //00000001
													   0x02,  //00000010
													   0x04,  //00000100
													   0x08,  //00001000
													   0x10,  //00010000
													   0x20,  //00100000
													   0x40,  //01000000
													   0x80   //10000000
};

//布隆过滤器参数
class bloom_parameters
{
public:
	//初始化参数
	bloom_parameters()
		: minimum_size(1),
		//std::numeric_limits模板库，存储数据类型的存储边界
		maximum_size(std::numeric_limits<unsigned long long int>::max()),
		minimum_number_of_hashes(1),
		maximum_number_of_hashes(std::numeric_limits<unsigned int>::max()),
		projected_element_count(10000),
		false_positive_probability(1.0 / projected_element_count),
		//ULL:unsigned long long
		random_seed(0xA5A5A5A55A5A5A5AULL)//1010 0101 ... 0101 1010
										  //  A    5        5    A
	{}

	virtual ~bloom_parameters()
	{}

	//重载!运算符，参数初始化不符合要求
	inline bool operator!()
	{
		return (minimum_size > maximum_size) ||
			(minimum_number_of_hashes > maximum_number_of_hashes) ||
			(minimum_number_of_hashes < 1) ||
			(0 == maximum_number_of_hashes) ||
			(0 == projected_element_count) ||
			(false_positive_probability < 0.0) ||
			//概率绝对值与double类型的正无穷大相等
			(std::numeric_limits<double>::infinity() == std::abs(false_positive_probability)) ||
			(0 == random_seed) ||
			(0xFFFFFFFFFFFFFFFFULL == random_seed);
	}

	// Allowable min/max size of the bloom filter in bits
	unsigned long long int minimum_size;
	unsigned long long int maximum_size;

	// Allowable min/max number of hash functions
	unsigned int minimum_number_of_hashes;
	unsigned int maximum_number_of_hashes;

	// The approximate number of elements to be inserted
	// into the bloom filter, should be within one order
	// of magnitude. The default is 10000.
	unsigned long long int projected_element_count;

	// The approximate false positive probability expected
	// from the bloom filter. The default is assumed to be
	// the reciprocal of the projected_element_count.
	double false_positive_probability;

	unsigned long long int random_seed;

	struct optimal_parameters_t
	{
		optimal_parameters_t()
			: number_of_hashes(0),
			table_size(0)
		{}

		unsigned int number_of_hashes;
		unsigned long long int table_size;
	};

	optimal_parameters_t optimal_parameters;

	virtual bool compute_optimal_parameters()
	{
		/*
		  Note:
		  The following will attempt to find the number of hash functions
		  and minimum amount of storage bits required to construct a bloom
		  filter consistent with the user defined false positive probability
		  and estimated element insertion count.
		*/

		//参数不符合要求，直接返回false
		if (!(*this))
			return false;

		double min_m = std::numeric_limits<double>::infinity();
		double min_k = 0.0;
		double k = 1.0;

		//误判率P=pow((1-exp(-k*n/m)),k),现在给定
		//误判率P,数据个数n，遍历哈希函数个数k，选择最小的位数组长度m
		//即m=-k*n/log(1-pow(p,1/k))
		while (k < 1000.0)
		{
			const double numerator = (-k * projected_element_count);
			const double denominator = std::log(1.0 - std::pow(false_positive_probability, 1.0 / k));

			const double curr_m = numerator / denominator;

			if (curr_m < min_m)
			{
				min_m = curr_m;
				min_k = k;
			}

			k += 1.0;
		}

		optimal_parameters_t& optp = optimal_parameters;

		optp.number_of_hashes = static_cast<unsigned int>(min_k);

		optp.table_size = static_cast<unsigned long long int>(min_m);
		//补足位数组长度为bits_per_char的倍数
		optp.table_size += (((optp.table_size % bits_per_char) != 0) ? (bits_per_char - (optp.table_size % bits_per_char)) : 0);

		//根据用户初始化参数选择hash函数个数和位数组长度，原则：小小取小（小于最小值取最小值），大大取大
		if (optp.number_of_hashes < minimum_number_of_hashes)
			optp.number_of_hashes = minimum_number_of_hashes;
		else if (optp.number_of_hashes > maximum_number_of_hashes)
			optp.number_of_hashes = maximum_number_of_hashes;

		if (optp.table_size < minimum_size)
			optp.table_size = minimum_size;
		else if (optp.table_size > maximum_size)
			optp.table_size = maximum_size;

		return true;
	}

};

//布隆过滤器
class bloom_filter
{
protected:

	typedef unsigned int bloom_type;
	typedef unsigned char cell_type;
	typedef std::vector<unsigned char> table_type;

public:

	bloom_filter()
		: salt_count_(0),
		table_size_(0),
		projected_element_count_(0),
		inserted_element_count_(0),
		random_seed_(0),
		desired_false_positive_probability_(0.0)
	{}

	bloom_filter(const bloom_parameters& p)
		: projected_element_count_(p.projected_element_count),
		inserted_element_count_(0),
		//random_seed(0xA5A5A5A55A5A5A5AULL)--question
		random_seed_((p.random_seed * 0xA5A5A5A5) + 1),
		desired_false_positive_probability_(p.false_positive_probability)
	{
		salt_count_ = p.optimal_parameters.number_of_hashes;
		table_size_ = p.optimal_parameters.table_size;

		generate_unique_salt();

		//重设位数组长度
		bit_table_.resize(table_size_ / bits_per_char, static_cast<unsigned char>(0x00));
	}

	bloom_filter(const bloom_filter& filter)
	{
		this->operator=(filter);
	}

	//重载相等判定符，判定两个filter参数是否一样
	inline bool operator == (const bloom_filter& f) const
	{
		if (this != &f)
		{
			return
				(salt_count_ == f.salt_count_) &&
				(table_size_ == f.table_size_) &&
				(bit_table_.size() == f.bit_table_.size()) &&
				(projected_element_count_ == f.projected_element_count_) &&
				(inserted_element_count_ == f.inserted_element_count_) &&
				(random_seed_ == f.random_seed_) &&
				(desired_false_positive_probability_ == f.desired_false_positive_probability_) &&
				(salt_ == f.salt_) &&
				(bit_table_ == f.bit_table_);
		}
		else
			return true;
	}

	//重载不等判定符，判定两个filter参数是否一样
	inline bool operator != (const bloom_filter& f) const
	{
		return !operator==(f);
	}

	//重载赋值符号，能够从另一个过滤器中复制相应的参数
	inline bloom_filter& operator = (const bloom_filter& f)
	{
		if (this != &f)
		{
			salt_count_ = f.salt_count_;
			table_size_ = f.table_size_;
			bit_table_ = f.bit_table_;
			salt_ = f.salt_;

			projected_element_count_ = f.projected_element_count_;
			inserted_element_count_ = f.inserted_element_count_;

			random_seed_ = f.random_seed_;

			desired_false_positive_probability_ = f.desired_false_positive_probability_;
		}

		return *this;
	}

	virtual ~bloom_filter()
	{}

	//返回元素个数是否为0，即位数组是否初始化
	inline bool operator!() const
	{
		return (0 == table_size_);
	}

	//清空位数组，设置已插入元素个数为0
	inline void clear()
	{
		std::fill(bit_table_.begin(), bit_table_.end(), static_cast<unsigned char>(0x00));
		inserted_element_count_ = 0;
	}

	//插入元素，提供const unsigned char*起始指针和长度
	inline void insert(const unsigned char* key_begin, const std::size_t& length)
	{
		std::size_t bit_index = 0;
		std::size_t bit = 0;

		//遍历随机数组构造哈希函数进行计算
		for (std::size_t i = 0; i < salt_.size(); ++i)
		{
			//计算哈希值对应的bit_index(第xx个数)，bit(位数组某一项的第x位)
			compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);
			//将位数组的第bit_index / bits_per_char位置的第bit位置1
			bit_table_[bit_index / bits_per_char] |= bit_mask[bit];
		}
		//插入元素个数加1
		++inserted_element_count_;
	}

	//插入模板，接受模板元素起始位置和长度
	template <typename T>
	inline void insert(const T& t)
	{
		// Note: T must be a C++ POD type.
		insert(reinterpret_cast<const unsigned char*>(&t), sizeof(T));
	}

	//插入元素，提供string
	inline void insert(const std::string& key)
	{
		insert(reinterpret_cast<const unsigned char*>(key.data()), key.size());
	}

	//插入元素，提供const char*
	inline void insert(const char* data, const std::size_t& length)
	{
		insert(reinterpret_cast<const unsigned char*>(data), length);
	}

	//输入模板，接受迭代器插入元素
	template <typename InputIterator>
	inline void insert(const InputIterator begin, const InputIterator end)
	{
		InputIterator itr = begin;

		while (end != itr)
		{
			insert(*(itr++));
		}
	}

	//检查元素是否在布隆过滤器中，接受元素起始位置和长度
	inline virtual bool contains(const unsigned char* key_begin, const std::size_t length) const
	{
		std::size_t bit_index = 0;
		std::size_t bit = 0;

		for (std::size_t i = 0; i < salt_.size(); ++i)
		{
			compute_indices(hash_ap(key_begin, length, salt_[i]), bit_index, bit);

			//只要有一个哈希函数计算结果不为1，返回false
			if ((bit_table_[bit_index / bits_per_char] & bit_mask[bit]) != bit_mask[bit])
			{
				return false;
			}
		}

		return true;
	}

	//检查模板，接受模板元素起始位置和长度
	template <typename T>
	inline bool contains(const T& t) const
	{
		return contains(reinterpret_cast<const unsigned char*>(&t), static_cast<std::size_t>(sizeof(T)));
	}

	//检查元素，提供string
	inline bool contains(const std::string& key) const
	{
		return contains(reinterpret_cast<const unsigned char*>(key.c_str()), key.size());
	}

	//检查元素，提供const char*
	inline bool contains(const char* data, const std::size_t& length) const
	{
		return contains(reinterpret_cast<const unsigned char*>(data), length);
	}

	//检查模板，接受迭代器检查全部的元素，
	//返回第一个不在过滤器中的元素迭代器对象或迭代器末尾对象（全在过滤器中）
	template <typename InputIterator>
	inline InputIterator contains_all(const InputIterator begin, const InputIterator end) const
	{
		InputIterator itr = begin;

		while (end != itr)
		{
			if (!contains(*itr))
			{
				return itr;
			}

			++itr;
		}

		return end;
	}

	//检查模板，接受迭代器检查全部的元素，
	//返回第一个在过滤器中的元素迭代器对象或迭代器末尾对象（全不在过滤器中）
	template <typename InputIterator>
	inline InputIterator contains_none(const InputIterator begin, const InputIterator end) const
	{
		InputIterator itr = begin;

		while (end != itr)
		{
			if (contains(*itr))
			{
				return itr;
			}

			++itr;
		}

		return end;
	}

	//返回位数组长度
	inline virtual unsigned long long int size() const
	{
		return table_size_;
	}

	//返回已插入元素个数
	inline unsigned long long int element_count() const
	{
		return inserted_element_count_;
	}

	//计算实际的误判率，用实际已插入的元素个数代替预设的元素个数
	inline double effective_fpp() const
	{
		/*
		  Note:
		  The effective false positive probability is calculated using the
		  designated table size and hash function count in conjunction with
		  the current number of inserted elements - not the user defined
		  predicated/expected number of inserted elements.
		*/
		//P=pow(1-exp(pow(-k*n/m)),k)
		return std::pow(1.0 - std::exp(-1.0 * salt_.size() * inserted_element_count_ / size()), 1.0 * salt_.size());
	}

	//重载&=符号，求两个过滤器的交集，即求两个过滤器共有的元素
	//前提：两个过滤器的哈希函数个数，位数组长度，随机种子相同
	inline bloom_filter& operator &= (const bloom_filter& f)
	{
		/* intersection */
		if (
			(salt_count_ == f.salt_count_) &&
			(table_size_ == f.table_size_) &&
			(random_seed_ == f.random_seed_)
			)
		{
			for (std::size_t i = 0; i < bit_table_.size(); ++i)
			{
				bit_table_[i] &= f.bit_table_[i];
			}
		}

		return *this;
	}

	//重载|=符号，求两个过滤器的并集，即合并两个过滤器中的元素
	//前提：两个过滤器的哈希函数个数，位数组长度，随机种子相同
	inline bloom_filter& operator |= (const bloom_filter& f)
	{
		/* union */
		if (
			(salt_count_ == f.salt_count_) &&
			(table_size_ == f.table_size_) &&
			(random_seed_ == f.random_seed_)
			)
		{
			for (std::size_t i = 0; i < bit_table_.size(); ++i)
			{
				bit_table_[i] |= f.bit_table_[i];
			}
		}

		return *this;
	}

	//重载^=符号，求两个过滤器的差集，即this过滤器除去f中的元素
	//前提：两个过滤器的哈希函数个数，位数组长度，随机种子相同
	inline bloom_filter& operator ^= (const bloom_filter& f)
	{
		/* difference */
		if (
			(salt_count_ == f.salt_count_) &&
			(table_size_ == f.table_size_) &&
			(random_seed_ == f.random_seed_)
			)
		{
			for (std::size_t i = 0; i < bit_table_.size(); ++i)
			{
				bit_table_[i] ^= f.bit_table_[i];
			}
		}

		return *this;
	}

	//返回位数组头指针
	inline const cell_type* table() const
	{
		return bit_table_.data();
	}

	//返回哈希函数个数
	inline std::size_t hash_count()
	{
		return salt_.size();
	}

protected:

	//根据哈希值计算需要置1的位数组index（即转化为第xx个数）和相应的bit位
	inline virtual void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const
	{
		bit_index = hash % table_size_;
		bit = bit_index % bits_per_char;
	}

	void generate_unique_salt()
	{
		/*
		  Note:
		  A distinct hash function need not be implementation-wise
		  distinct. In the current implementation "seeding" a common
		  hash function with different values seems to be adequate.
		*/
		const unsigned int predef_salt_count = 128;

		static const bloom_type predef_salt[predef_salt_count] =
		{
		   0xAAAAAAAA, 0x55555555, 0x33333333, 0xCCCCCCCC,
		   0x66666666, 0x99999999, 0xB5B5B5B5, 0x4B4B4B4B,
		   0xAA55AA55, 0x55335533, 0x33CC33CC, 0xCC66CC66,
		   0x66996699, 0x99B599B5, 0xB54BB54B, 0x4BAA4BAA,
		   0xAA33AA33, 0x55CC55CC, 0x33663366, 0xCC99CC99,
		   0x66B566B5, 0x994B994B, 0xB5AAB5AA, 0xAAAAAA33,
		   0x555555CC, 0x33333366, 0xCCCCCC99, 0x666666B5,
		   0x9999994B, 0xB5B5B5AA, 0xFFFFFFFF, 0xFFFF0000,
		   0xB823D5EB, 0xC1191CDF, 0xF623AEB3, 0xDB58499F,
		   0xC8D42E70, 0xB173F616, 0xA91A5967, 0xDA427D63,
		   0xB1E8A2EA, 0xF6C0D155, 0x4909FEA3, 0xA68CC6A7,
		   0xC395E782, 0xA26057EB, 0x0CD5DA28, 0x467C5492,
		   0xF15E6982, 0x61C6FAD3, 0x9615E352, 0x6E9E355A,
		   0x689B563E, 0x0C9831A8, 0x6753C18B, 0xA622689B,
		   0x8CA63C47, 0x42CC2884, 0x8E89919B, 0x6EDBD7D3,
		   0x15B6796C, 0x1D6FDFE4, 0x63FF9092, 0xE7401432,
		   0xEFFE9412, 0xAEAEDF79, 0x9F245A31, 0x83C136FC,
		   0xC3DA4A8C, 0xA5112C8C, 0x5271F491, 0x9A948DAB,
		   0xCEE59A8D, 0xB5F525AB, 0x59D13217, 0x24E7C331,
		   0x697C2103, 0x84B0A460, 0x86156DA9, 0xAEF2AC68,
		   0x23243DA5, 0x3F649643, 0x5FA495A8, 0x67710DF8,
		   0x9A6C499E, 0xDCFB0227, 0x46A43433, 0x1832B07A,
		   0xC46AFF3C, 0xB9C8FFF0, 0xC9500467, 0x34431BDF,
		   0xB652432B, 0xE367F12B, 0x427F4C1B, 0x224C006E,
		   0x2E7E5A89, 0x96F99AA5, 0x0BEB452A, 0x2FD87C39,
		   0x74B2E1FB, 0x222EFD24, 0xF357F60C, 0x440FCB1E,
		   0x8BBE030F, 0x6704DC29, 0x1144D12F, 0x948B1355,
		   0x6D8FD7E9, 0x1C11A014, 0xADD1592F, 0xFB3C712E,
		   0xFC77642F, 0xF9C4CE8C, 0x31312FB9, 0x08B0DD79,
		   0x318FA6E7, 0xC040D23D, 0xC0589AA7, 0x0CA5C075,
		   0xF874B172, 0x0CF914D5, 0x784D3280, 0x4E8CFEBC,
		   0xC569F575, 0xCDB2A091, 0x2CC016B4, 0x5C5F4421
		};

		if (salt_count_ <= predef_salt_count)
		{
			//左闭右开，第三个参数是提供另一个有序组的头地址
			std::copy(predef_salt,
				predef_salt + salt_count_,
				//插入型迭代器
				std::back_inserter(salt_));

			for (std::size_t i = 0; i < salt_.size(); ++i)
			{
				/*
				   Note:
				   This is done to integrate the user defined random seed,
				   so as to allow for the generation of unique bloom filter
				   instances.
				*/
				//question
				salt_[i] = salt_[i] * salt_[(i + 3) % salt_.size()] + static_cast<bloom_type>(random_seed_);
			}
		}
		else
		{
			//程序中自定义的salt(128个)不够
			std::copy(predef_salt, predef_salt + predef_salt_count, std::back_inserter(salt_));

			srand(static_cast<unsigned int>(random_seed_));

			while (salt_.size() < salt_count_)
			{
				bloom_type current_salt = static_cast<bloom_type>(rand()) * static_cast<bloom_type>(rand());

				if (0 == current_salt)
					continue;

				//salt_中最后一个数与当前生成的不同时，插入，即没有位置相邻数值相同的随机数
				if (salt_.end() == std::find(salt_.begin(), salt_.end(), current_salt))
				{
					salt_.push_back(current_salt);
				}
			}
		}
	}

	//哈希值计算，提供元素起始地址，长度和哈希类型（随机数组中的一项）
	inline bloom_type hash_ap(const unsigned char* begin, std::size_t remaining_length, bloom_type hash) const
	{
		const unsigned char* itr = begin;
		unsigned int loop = 0;

		//元素长度大于8个字节
		while (remaining_length >= 8)
		{
			//unsigned int 4个字节，获取元素的前八个字节并分解为两个unsigned int
			const unsigned int& i1 = *(reinterpret_cast<const unsigned int*>(itr)); itr += sizeof(unsigned int);
			const unsigned int& i2 = *(reinterpret_cast<const unsigned int*>(itr)); itr += sizeof(unsigned int);

			//四字节（32位）的哈希数，question:哈希函数定义
			hash ^= (hash << 7) ^ i1 * (hash >> 3) ^
				(~((hash << 11) + (i2 ^ (hash >> 5))));

			remaining_length -= 8;
		}

		//元素长度小于8或已有部分进行hash,剩余部分长度小于8
		if (remaining_length)
		{
			//元素长度大于4，截取前四个元素作为unsigned int,根据是前四个还是后四个选择相应的哈希函数
			if (remaining_length >= 4)
			{
				const unsigned int& i = *(reinterpret_cast<const unsigned int*>(itr));

				if (loop & 0x01)
					hash ^= (hash << 7) ^ i * (hash >> 3);
				else
					hash ^= (~((hash << 11) + (i ^ (hash >> 5))));

				++loop;

				remaining_length -= 4;

				itr += sizeof(unsigned int);
			}
			//元素长度大于2，截取前2个元素作为unsigned short,根据是前四个还是后四个选择相应的哈希函数
			if (remaining_length >= 2)
			{
				const unsigned short& i = *(reinterpret_cast<const unsigned short*>(itr));

				if (loop & 0x01)
					hash ^= (hash << 7) ^ i * (hash >> 3);
				else
					hash ^= (~((hash << 11) + (i ^ (hash >> 5))));

				++loop;

				remaining_length -= 2;

				itr += sizeof(unsigned short);
			}
			//元素长度为1
			if (remaining_length)
			{
				hash += ((*itr) ^ (hash * 0xA5A5A5A5)) + loop;
			}
		}

		return hash;
	}

	std::vector<bloom_type>    salt_;	//随机数组（哈希函数相关）
	std::vector<unsigned char> bit_table_;	//位数组
	unsigned int               salt_count_;	//随机数组长度
	unsigned long long int     table_size_;	//元素个数
	unsigned long long int     projected_element_count_;	//预设元素个数
	unsigned long long int     inserted_element_count_;		//已插入元素个数
	unsigned long long int     random_seed_;	//随机种子
	double                     desired_false_positive_probability_;	//误判率
};

//重载&运算符，返回两个布隆过滤器的交集
inline bloom_filter operator & (const bloom_filter& a, const bloom_filter& b)
{
	bloom_filter result = a;
	result &= b;
	return result;
}

//重载|运算符，返回两个布隆过滤器的并集
inline bloom_filter operator | (const bloom_filter& a, const bloom_filter& b)
{
	bloom_filter result = a;
	result |= b;
	return result;
}

//重载^运算符，返回两个布隆过滤器的差集
inline bloom_filter operator ^ (const bloom_filter& a, const bloom_filter& b)
{
	bloom_filter result = a;
	result ^= b;
	return result;
}

//public继承布隆过滤器，得到可压缩的布隆过滤器
class compressible_bloom_filter : public bloom_filter
{
public:
	//通过布隆过滤器参数构造布隆过滤器，并将其位长度转入位长度列表
	compressible_bloom_filter(const bloom_parameters& p)
		: bloom_filter(p)
	{
		size_list.push_back(table_size_);
	}

	//获得最后一个位长度
	inline unsigned long long int size() const
	{
		return size_list.back();
	}


	inline bool compress(const double& percentage)
	{
		if (
			(percentage < 0.0) ||
			(percentage >= 100.0)
			)
		{
			return false;
		}

		//获得原始位长度，使用压缩率计算新的位长度
		unsigned long long int original_table_size = size_list.back();
		unsigned long long int new_table_size = static_cast<unsigned long long int>((size_list.back() * (1.0 - (percentage / 100.0))));

		//使位长度能被bits_per_char整除
		new_table_size -= new_table_size % bits_per_char;

		//新的位长度比bits_per_char还小或比原来位长度还大（无法压缩）
		if (
			(bits_per_char > new_table_size) ||
			(new_table_size >= original_table_size)
			)
		{
			return false;
		}

		//计算实际的误判率
		desired_false_positive_probability_ = effective_fpp();

		//使用压缩后的位长度构建新的位数组，并将原位数组的前面一部分拷贝进去
		const unsigned long long int new_tbl_raw_size = new_table_size / bits_per_char;

		table_type tmp(new_tbl_raw_size);

		std::copy(bit_table_.begin(), bit_table_.begin() + new_tbl_raw_size, tmp.begin());

		//使用迭代器拷贝原位数组的后面部分，与前面部分求并，question(后半部分较大，产生内存覆盖)
		typedef table_type::iterator itr_t;

		itr_t itr = bit_table_.begin() + (new_table_size / bits_per_char);
		itr_t end = bit_table_.begin() + (original_table_size / bits_per_char);
		itr_t itr_tmp = tmp.begin();

		while (end != itr)
		{
			*(itr_tmp++) |= (*itr++);
		}

		std::swap(bit_table_, tmp);

		size_list.push_back(new_table_size);

		return true;
	}

private:

	//重写父类中哈希值计算需要置1的位数组index（即转化为第xx个数）和相应的bit位
	inline void compute_indices(const bloom_type& hash, std::size_t& bit_index, std::size_t& bit) const
	{
		bit_index = hash;

		//对位长度列表中的每个位长度循环求余，可得到该哈希值的相应位置
		for (std::size_t i = 0; i < size_list.size(); ++i)
		{
			bit_index %= size_list[i];
		}

		bit = bit_index % bits_per_char;
	}

	std::vector<unsigned long long int> size_list;	//位长度存储列表
};

#endif


/*
  Note 1:
  If it can be guaranteed that bits_per_char will be of the form 2^n then
  the following optimization can be used:

  bit_table_[bit_index >> n] |= bit_mask[bit_index & (bits_per_char - 1)];

  Note 2:
  For performance reasons where possible when allocating memory it should
  be aligned (aligned_alloc) according to the architecture being used.
*/
