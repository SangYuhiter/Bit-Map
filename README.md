# Bit-Map-and-Bloom-Filter
位图算法和布隆过滤器的实现

## 实验内容1：在100万个不重复的正整数中，是否存在某个数。分别使用位图和bloomfilter实现。

### 1 ：生成1万个随机数 （srand  rand函数 随机数范围0到rand_max）， 
### 2 ：用bitmap表示这1万个随机数，
### 3：将1万个随机数放到bloomfilter容器中
### 4：再产生1000个随机数，通过bitmap和bloomfilter判断，检测bloomfilter的误判率

## 实验内容2:  通过给定的数据集，两个分别具有大量URL地址信息的数据文件，找出两个文件共同的URL。 

