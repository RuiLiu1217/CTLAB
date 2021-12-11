#ifndef PHANTOM_CREATOR_H
#define PHANTOM_CREATOR_H
#include <vector>
#include <string>
class PhantomCreator {
public:
	PhantomCreator();
	void create(int sizex, int sizey, int sizez);
	void reset();
	bool save(std::string& fileName);
private:
	std::vector<float> data;
};

#endif