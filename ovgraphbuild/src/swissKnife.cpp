#include "swissKnife.h"

template<typename T>
std::vector<T> intersection(std::vector<T> &v1, std::vector<T> &v2)
{
    std::vector<T> v3;

    std::sort(std::begin(v1), std::end(v1));
    std::sort(std::begin(v2), std::end(v2));

    std::set_intersection(std::begin(v1),
                          std::end(v1),
                          std::begin(v2),
                          std::end(v2),
                          std::back_inserter(v3));

    return v3;
}
