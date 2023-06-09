#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>     
#include <chrono> 
#include <random>
#include <algorithm>
#include "Tools.cpp"

#define MUTATION_PROB 0.2

using namespace std::chrono;
using namespace std;


class NumberCity{
  public:
     
       int index = 0;
      float x = 0.0;
      float y = 0.0;

      NumberCity(int index, float xcoord, float ycoord){this->index = index;this->x=xcoord; this->y=ycoord;};
      NumberCity() = default;
      ~NumberCity(){x=0.0;y=0.0;};

      float calcDistance(NumberCity next){
          
          float xStart = abs(this->x - next.x);
          float yStart = abs(this->y - next.y);

          return sqrt(pow(xStart, 2) + pow(yStart, 2));
      }
};


class Cromozom    // si vom aveam o populatie de astfel de cromozomi
{
    public:
      vector<NumberCity> route; 
      size_t route_dimension;       //route dimension (how many cities are)
      float fitness = 0.0;
      float distance = 0.0;
      float rprob =0.0;       //prob cromozomului in cross-over

      Cromozom(size_t dim){route_dimension = dim;};
      Cromozom() = default;
      ~Cromozom(){route.clear();};

      float getFitness()
      {
          return 1 / this->getDistance();
      }

      float getDistance(){
            this->distance = 0.0;
            for(int i = 1; i < route_dimension; ++i)
                distance += route[i].calcDistance(route[i-1]);

            return this->distance;
      }

      void random_route()
      {
        /* shuffles the city order inside chromozome */
        random_device random;
        std::shuffle(route.begin(), route.end(), random);
      }

     void RandomProb()
     {
         rprob = generateFloatRandom(0, 1);
     } 

     void mutate()
     {
        int index = 0;
        for(int i = 0; i < this->route_dimension; i++)
        {
            float r = generateFloatRandom(0, 1);
            if(r < MUTATION_PROB)
            {
                index = generateRandom(0, this->route_dimension);
                while(index == i)
                {
                    index = generateRandom(0, this->route_dimension);
                }
                swap(this->route[i], this->route[index]);
            }
        }
     }

     void update_fitness()
     {
        this->fitness = this->getFitness();
     }

     void swapp_hemming()
     {
         int index1 = generateRandom(0, this->route_dimension);
         int index2 = generateRandom(0, this->route_dimension);
         while(index1 == index2)
         {
            index2 = generateRandom(0, this->route_dimension);
         }
         swap(this->route[index1], this->route[index2]);
     }

     void rgibnm()
     {  
        int index1=generateRandom(0,this->route_dimension);
        int index2=0;
        float min=1000000;
        float k;
       //Cautam orasul de distanta minima
        for(int i=0;i<this->route_dimension;i++)
        {
           if(i!=index1)
           {
              k = this->route[i].calcDistance(this->route[index1]);  //distanta intre orasul i si orasul index 1;
               if(min>k)//Distanta dintre 2 orase calculez
               {
                  min=k;
                  index2=i;
               }
           }
        }
       // cout<<k<<'\n';
        //normal trebuie sa caut alt oras apropiat de orasul cu index ul 2;
        
        if(index1 < this->route_dimension-1)
                index1++;
                else
                index1=0;

        if(index1 != index2)
        { swap(this->route[index1],this->route[index2]);
         //cout<<index1<<'\n';
         //cout<<index2<<'\n';
        }
     }
      void swapp_reverse()
      {
         int index1 = generateRandom(0, this->route_dimension);
         int index2 = generateRandom(0, this->route_dimension); 
        
        while(index1 == index2)
         {
            index2 = generateRandom(0, this->route_dimension);
         }

         reverse(this->route.begin() + index1, this->route.begin() + index2);
      }

     void apply_SA()
     {
        Cromozom neighbor;
        Cromozom tmp;
        tmp.route = this->route;
        tmp.route_dimension = this->route_dimension;
        float Cost1=0;  
        float Cost2=0;
        Cost1=tmp.getDistance();     
        float t=100;
        while(t>0.0000001)
        {
            neighbor=tmp;
            //neighbor.swapp_hemming(); 
            neighbor.swapp_reverse();
           // neighbor.rgibnm();
            //neighbor.mutate();                    
            Cost2=neighbor.getDistance(); 
            if(Cost2<Cost1)
            {
                tmp=neighbor;
                Cost1=Cost2;
            }
            else
            {
                float abso = abs(Cost2 - Cost1);
                float exponent = exp( -abso / t );
                float random = generateFloatRandom(0, 1 - 0.000001);
                if(random < exponent)
                {
                    tmp = neighbor;
                    Cost1=Cost2;
                    t=t*0.9995;
                }
            }
        }
        this->route = tmp.route;


        cout << "distance: " << this->getDistance() << endl;
    }

};



bool orderByFitness(const Cromozom &a, const Cromozom &b) {
        return a.fitness > b.fitness;
}
 
bool orderByProb(const Cromozom &a, const Cromozom &b){
       return a.rprob < b.rprob;
}

Cromozom Read(char* source){
    ifstream in;
    in.open(source);
    size_t dim = 0;
    in>>dim;
    string line;
    float coordinate_x = 0.0, coordinate_y = 0.0;
    int index = 0;
    Cromozom initial_crom(dim);
    while(in >> index){
        in >> coordinate_x;
        in >> coordinate_y;
        NumberCity c{index, coordinate_x, coordinate_y};
        initial_crom.route.emplace_back(c);
    }
    return initial_crom;
}

void Read_best(char* source, vector<int> &order){
    ifstream in;
    in.open(source);
    int index = 0;
    while(in >> index){
        order.emplace_back(index);
    }
}