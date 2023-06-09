#include <iostream>
#include <chrono>
#include <random>

using namespace std;
using namespace std::chrono;

mt19937 generator(0);

long long get_current_time() {
    milliseconds ms = duration_cast< milliseconds >(
        system_clock::now().time_since_epoch()
    );
    return ms.count();
}

float generateFloatRandom(float x,float y)
{
     uniform_real_distribution<float> distribution(x, y); // [x, y]
    long long ct = get_current_time() % 10000;

    float rn = distribution(generator) * ct;
    float dif = y - x;

    while(rn >= dif) 
    {
        rn -= dif;
    }

    rn += x;

    return rn;
}

int generateRandom(int x, int y) {
    uniform_int_distribution<int> distribution(x, y); // [x, y]
    long long ct = get_current_time() % 10000;
    int rn = distribution(generator) * ct;
    
    int dif = y - x;
    return (rn % dif) + x;
}

float Ma(float vec[],int n)
{
    float S = 0;
    
    for (int i = 0; i < n; i++)
        S = S + vec[i];
    S = S / n;
    return S ;
   
}

float standardDeviation(float vect[],int n)
{
    float Sum = 0;
    float a = Ma(vect,n);
    for (int i = 0; i < n; i++)
    {
        float  x = vect[i];
        Sum = Sum + (x - a)*(x - a);
    }
    return sqrt(Sum / n);
}
