#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <unordered_map>

#include"Cromozom.cpp"

#define DimPop 200
#define M_C_PROB 0.45

using namespace std;


void Mutate(vector<Cromozom> &population)   //am un vect cu traseul actual.
{
  for(int i = 2; i < DimPop; i++)
  {
   // population[i].mutate();
   // population[i].swapp_hemming();
    float r = generateFloatRandom(0, 1);
      if(r < MUTATION_PROB)
        {
             population[i].swapp_hemming();
       }
  }
}

void Inverse(vector<Cromozom> &population)
{
 for(int i = 2; i < DimPop; i++)
  {
   // population[i].mutate();
   // population[i].swapp_hemming();
    float r = generateFloatRandom(0, 1);
      if(r < MUTATION_PROB)
        {
             population[i].swapp_reverse();
       }
  }
}

/*
void HyperMutation(vector<Cromozom> &population)   //am un vect cu traseul actual.
{
  for(int i = 2; i < DimPop; i++)
  {
   // population[i].mutate();
   // population[i].swapp_hemming();
     float r = generateFloatRandom(0, 1);
      if(r < MUTATION_PROB*2.5)
        {
             population[i].swapp_reverse();
       }
  }
}
*/
void RGIBNM(vector<Cromozom> &population)
{
      for(int i = 2; i < DimPop; i++)
  {
   // population[i].mutate();
   // population[i].swapp_hemming();
    float r = generateFloatRandom(0, 1);
      if(r < MUTATION_PROB)
        {
             population[i].rgibnm();
       }
  }
}

void IRGIBNM(vector<Cromozom> &population) //Improved
{
     for(int i = 2; i < DimPop; i++)
  {
    float r = generateFloatRandom(0, 1);
      if(r < MUTATION_PROB)
        {
             population[i].swapp_reverse();
             population[i].rgibnm();
       }
  }
}

void SBM(vector<Cromozom> &population)   //Select Best Mutation
{
     vector<Cromozom> &populationClone1 = population;
     vector<Cromozom> &populationClone2 = population;
     vector<Cromozom> &populationClone3 = population;
     vector<Cromozom> &populationClone4 = population;
    
    for(int i = 2; i < DimPop; i++)
     { 
      float r = generateFloatRandom(0, 1);
       if(r < MUTATION_PROB)
        {
            float a = population[i].getDistance();

             populationClone1[i].swapp_hemming();
              float b=populationClone1[i].getDistance();
            
             populationClone2[i].swapp_reverse();
              float c=populationClone2[i].getDistance();
            
             populationClone3[i].rgibnm();
              float d=populationClone3[i].getDistance(); 
            
             populationClone4[i].swapp_reverse();
             populationClone4[i].rgibnm();
               float e=populationClone4[i].getDistance();
            
             //cout<<a<<' '<<b<<' '<<c<<' '<<d<<' '<<e<<'\n';
            
             if(a < b)
                {
                  a=b;
                  population[i]=populationClone1[i];
                  //cout<<"da"<<'\n';
                }
             
             if(a < c)
               {
                  a=c;
                  population[i]=populationClone2[i]; 
                 // cout<<"da"<<'\n';
               }
             
             if(a < d)
              {
                  a=d;
                  population[i]=populationClone3[i];
                //  cout<<"da"<<'\n';
              }
             
             if(a < e)
             {
                  a=e;
                  population[i]=populationClone4[i];
                  //cout<<"da"<<'\n';
             }
      
       }
     } 
}

void ElitOXCross(Cromozom &d1, Cromozom &d2,  Cromozom &p1, Cromozom &p2)
{
      unordered_map<int, NumberCity>order_p1, order_p2;
    NumberCity c{}, c1{};
    int index_1=0, index_2=0;
    int selected_1[p1.route_dimension+2], selected_2[p2.route_dimension+2];


    for(int i = 0; i < p1.route_dimension - 1; i++)
    {
        // we create the order by knowing what is the next city from the current city 
        order_p1.insert(make_pair(p1.route[i].index, p1.route[i+1]));
        order_p2.insert(make_pair(p2.route[i].index, p2.route[i+1]));
        selected_1[i] = 0;
        selected_2[i] = 0;
    }
    order_p1.insert(make_pair(p1.route[p1.route_dimension-1].index, p1.route[0]));
    order_p2.insert(make_pair(p2.route[p2.route_dimension-1].index, p2.route[0]));

    selected_1[p1.route_dimension-1] = selected_2[p1.route_dimension-1] = 0;
    selected_1[p1.route_dimension] = selected_2[p1.route_dimension] = 0;

    // create d1 
    d1.route.emplace_back(p2.route[0]);
    selected_1[p2.route[0].index] = 1;
     // create d2 
    d2.route.emplace_back(p1.route[0]);
    selected_2[p1.route[0].index] = 1;
    for(int i = 1; i <= p1.route_dimension - 1; i++)
    {
        c = order_p1[d1.route[index_1].index];
        c1 = order_p2[d1.route[index_1].index];
        float distance_between_d1_and_closest_city_p1 = d1.route[index_1].calcDistance(c);
        float distance_between_d1_and_closest_city_p2 = d1.route[index_1].calcDistance(c1);
        if(
            (distance_between_d1_and_closest_city_p2 < distance_between_d1_and_closest_city_p1) &&
            (selected_1[c.index] == 0)
        )
        {
            d1.route.emplace_back(c); 
            selected_1[c.index] = 1;
            index_1++;
        }
        else if(selected_1[c1.index] == 0)
        {
            d1.route.emplace_back(c1);
            selected_1[c1.index] = 1;
            index_1++;
        }
        else
        {
            for(auto const& val : order_p2)
            {
                {
                    if(selected_1[val.second.index] == 0)
                    {                    
                        d1.route.emplace_back(val.second);
                        selected_1[val.second.index] = 1;
                        index_1++;
                        break;
                    }

                }  
            }                   
        }

        c1 = order_p2[d2.route[index_2].index];
        c = order_p1[d2.route[index_2].index];
        float distance_between_d2_and_closest_city_p1 = d2.route[index_2].calcDistance(c);
        float distance_between_d2_and_closest_city_p2 = d2.route[index_2].calcDistance(c1);
        if(
            (distance_between_d2_and_closest_city_p1 < distance_between_d2_and_closest_city_p2) &&
            (selected_2[c.index] == 0)
        )
        {
            d2.route.emplace_back(c); 
            selected_2[c.index] = 1;
            index_2++;
        }
        else if(selected_2[c1.index] == 0)
        {
            d2.route.emplace_back(c1);
            selected_2[c1.index] = 1;
            index_2++;
        }
        else
        {
            for(auto const& val : order_p1)
            {
                {
                    if(selected_2[val.second.index] == 0)
                    {                    
                        d2.route.emplace_back(val.second);
                        selected_2[val.second.index] = 1;
                        index_2++;
                    }

                }  
            }                   
        }

    }
}
void ElitistPopCross1(vector<Cromozom> &pop)
{
    Cromozom best_crom = pop[0];
    for (int i = 0; i < DimPop; i++) {
        pop[i].RandomProb();
    }
    sort(pop.begin(), pop.end(), orderByProb);
    for(int i = 0; i < 20; i++)
        pop[i] = best_crom;

    int cutting_point_for_pop = 0;
    while (pop[cutting_point_for_pop].rprob < M_C_PROB && cutting_point_for_pop < DimPop) {
        cutting_point_for_pop++;
    }
    cutting_point_for_pop--;

    if (cutting_point_for_pop % 2 == 1)
        cutting_point_for_pop += (generateRandom(1, 10) % 2 ? 1 : -1);


    for (int i = 2; i < cutting_point_for_pop; i+=2) {
        Cromozom d1{pop[i].route_dimension}, d2{pop[i+1].route_dimension};
        ElitOXCross(d1, d2, pop[i], pop[i+1]);
        pop[i] = d1;
        pop[i+1] = d2;
    }
}

void OXCross(Cromozom &d1,  Cromozom &p1, Cromozom &p2)   //asta merge cu 2 taieturi desi m as gadi sa l fac cu un nr random pana in m/2 taieturi (iar bucatiile dintre acele taietuti)                                               //raman la fel.
{
    vector<NumberCity>route(p1.route_dimension);
    unordered_map<int, NumberCity> visited;
    for(int i = 0; i < p1.route_dimension; i++)
        route[i] = NumberCity{};

    int first_c_point = (int)(generateFloatRandom(0, 1) * p1.route_dimension);
    int second_c_point = (int)(generateFloatRandom(0, 1) * p1.route_dimension);

    for(int i = 0; i < p1.route_dimension; i++)
    {
        if(first_c_point < second_c_point && i > first_c_point && i < second_c_point)
        {
            route[i]=p1.route[i];
            visited.insert(make_pair(p1.route[i].index, p1.route[i]));
        }
        else if(first_c_point > second_c_point)
        {
            if(!(i < first_c_point && i > second_c_point))
            {
                route[i]=p1.route[i];
                visited.insert(make_pair(p1.route[i].index, p1.route[i]));
            }
        }
    }

    for(int i = 0; i < p2.route_dimension; i++)
    {
        if(visited.find(p2.route[i].index) == visited.end())
        {
            for(int j = 0; j < p2.route_dimension;j++)
            {
                if(route[j].index == 0)
                {
                    route[j] = p2.route[i];
                    break;
                }
            }
        }
    }
    d1.route = route;
    d1.route_dimension = p1.route_dimension;
}

void OXCrossover(vector<Cromozom> &pop)
{
    Cromozom best_crom = pop[0];
    for (int i = 0; i < DimPop; i++) {
        pop[i].RandomProb();
    }
    sort(pop.begin(), pop.end(), orderByProb);
    for(int i = 0; i < 20; i++)
        pop[i] = best_crom;

    int cutting_point_for_pop = 0;
    while (pop[cutting_point_for_pop].rprob < M_C_PROB && cutting_point_for_pop < DimPop) {
        cutting_point_for_pop++;
    }
    cutting_point_for_pop--;

    if (cutting_point_for_pop % 2 == 1)
        cutting_point_for_pop += (generateRandom(1, 10) % 2 ? 1 : -1);
    
    
    for (int i = 2; i < cutting_point_for_pop; i++) 
    {
        Cromozom d1{pop[i].route_dimension};
        OXCross(d1, pop[i], pop[i+1]);
        pop[i] = d1;
    } 
}

void CXCross(Cromozom &d1, Cromozom &d2,  Cromozom &p1, Cromozom &p2)
{
  d1.route[0] = p1.route[0];
  d2.route[0] = p2.route[0]; 
  int a,b;
  for(int i=1;i<p1.route_dimension;i++)
   {
      a = p1.route[i-1].index;
      b = p2.route[a].index;
      d1.route[i]=p1.route[b];

      a = p2.route[i-1].index;
      b = p1.route [a].index;
      d2.route[i]=p2.route[b];
   }
}

void CXCrossOver(vector<Cromozom> &pop)
{
     Cromozom best_crom = pop[0];
    for (int i = 0; i < DimPop; i++) {
        pop[i].RandomProb();
    }
    sort(pop.begin(), pop.end(), orderByProb);
   /* for(int i = 0; i < 20; i++)
        pop[i] = best_crom;
   */
    int cutting_point_for_pop = 0;
    while (pop[cutting_point_for_pop].rprob < M_C_PROB && cutting_point_for_pop < DimPop) {
        cutting_point_for_pop++;
    }
    cutting_point_for_pop--;

    if (cutting_point_for_pop % 2 == 1)
        cutting_point_for_pop += (generateRandom(1, 10) % 2 ? 1 : -1);


    for (int i = 2; i < cutting_point_for_pop; i+=2) {
        Cromozom d1{pop[i].route_dimension}, d2{pop[i+1].route_dimension};
        CXCross(d1, d2, pop[i], pop[i+1]);
        pop[i] = d1;
        pop[i+1] = d2;
    }
}

void IGXCross(Cromozom &d1, Cromozom &d2,  Cromozom &p1, Cromozom &p2)
{
     unordered_map<int, NumberCity>order_p1, order_p2;
    NumberCity c{}, c1{};
    int index_1=0, index_2=0;
    int selected_1[p1.route_dimension+2], selected_2[p2.route_dimension+2];


    for(int i = 0; i < p1.route_dimension - 1; i++)
    {
        // we create the order by knowing what is the next city from the current city 
        order_p1.insert(make_pair(p1.route[i].index, p1.route[i+1]));
        order_p2.insert(make_pair(p2.route[i].index, p2.route[i+1]));
        selected_1[i] = 0;
        selected_2[i] = 0;
    }
    order_p1.insert(make_pair(p1.route[p1.route_dimension-1].index, p1.route[0]));
    order_p2.insert(make_pair(p2.route[p2.route_dimension-1].index, p2.route[0]));

    selected_1[p1.route_dimension-1] = selected_2[p1.route_dimension-1] = 0;
    selected_1[p1.route_dimension] = selected_2[p1.route_dimension] = 0;

    // create d1 
    d1.route.emplace_back(p2.route[0]);
    selected_1[p2.route[0].index] = 1;
     // create d2 
    d2.route.emplace_back(p1.route[0]);
    selected_2[p1.route[0].index] = 1;
    for(int i = 1; i <= p1.route_dimension - 1; i++)
    {
        c = order_p1[d1.route[index_1].index];
        c1 = order_p2[d1.route[index_1].index];
        float distance_between_d1_and_closest_city_p1 = d1.route[index_1].calcDistance(c);
        float distance_between_d1_and_closest_city_p2 = d1.route[index_1].calcDistance(c1);
        if(
            (distance_between_d1_and_closest_city_p2 < distance_between_d1_and_closest_city_p1) &&
            (selected_1[c.index] == 0)
        )
        {
            d1.route.emplace_back(c); 
            selected_1[c.index] = 1;
            index_1++;
        }
        else if(selected_1[c1.index] == 0)
        {
            d1.route.emplace_back(c1);
            selected_1[c1.index] = 1;
            index_1++;
        }
        else
        {
            for(auto const& val : order_p2)
            {
                {
                    if(selected_1[val.second.index] == 0)
                    {                    
                        d1.route.emplace_back(val.second);
                        selected_1[val.second.index] = 1;
                        index_1++;
                        break;
                    }

                }  
            }                   
        }

        c1 = order_p2[d2.route[index_2].index];
        c = order_p1[d2.route[index_2].index];
        float distance_between_d2_and_closest_city_p1 = d2.route[index_2].calcDistance(c);
        float distance_between_d2_and_closest_city_p2 = d2.route[index_2].calcDistance(c1);
        if(
            (distance_between_d2_and_closest_city_p1 < distance_between_d2_and_closest_city_p2) &&
            (selected_2[c.index] == 0)
        )
        {
            d2.route.emplace_back(c); 
            selected_2[c.index] = 1;
            index_2++;
        }
        else if(selected_2[c1.index] == 0)
        {
            d2.route.emplace_back(c1);
            selected_2[c1.index] = 1;
            index_2++;
        }
        else
        {
            for(auto const& val : order_p1)
            {
                {
                    if(selected_2[val.second.index] == 0)
                    {                    
                        d2.route.emplace_back(val.second);
                        selected_2[val.second.index] = 1;
                        index_2++;
                    }

                }  
            }                   
        }

    }
}
 


void IGXCrossOver(vector<Cromozom> &pop)
{
    Cromozom best_crom = pop[0];
    for (int i = 0; i < DimPop; i++) {
        pop[i].RandomProb();
    }
    sort(pop.begin(), pop.end(), orderByProb);
    for(int i = 0; i < 20; i++)
        pop[i] = best_crom;

    int cutting_point_for_pop = 0;
    while (pop[cutting_point_for_pop].rprob < M_C_PROB && cutting_point_for_pop < DimPop) {
        cutting_point_for_pop++;
    }
    cutting_point_for_pop--;

    if (cutting_point_for_pop % 2 == 1)
        cutting_point_for_pop += (generateRandom(1, 10) % 2 ? 1 : -1);


    for (int i = 2; i < cutting_point_for_pop; i+=2) {
        Cromozom d1{pop[i].route_dimension}, d2{pop[i+1].route_dimension};
        IGXCross(d1, d2, pop[i], pop[i+1]);
        pop[i] = d1;
        pop[i+1] = d2;
    }
}


void SelectWheelOfFortune(vector<Cromozom> &population)
{
  //AICI facem ruleta dar combinat cu un elitism (salvez numai una care e cea mai buna in toate gen)
  //Va trebui sa vedem cand incepe ca convearga cam tare si atunci la restul oraselor in afara de primul va trb sa le regenerez random ; ininte de fiecare selectie in functia geneticului verific sd.
    float SelectionWheel[DimPop * 2 + 1];  //de 401
    float TotalScore=0;
    SelectionWheel[0] = 0;

    for(int i=0;i<DimPop;i++)
        TotalScore += population[i].getFitness();

    for(int i=1;i<=DimPop;i++)
        SelectionWheel[i] = population[i - 1].fitness/TotalScore;

    for(int i=1;i<=DimPop;i++)
        SelectionWheel[i]+=SelectionWheel[i-1];

     for(int i = 2; i<DimPop; i++)
     {
        float r= generateFloatRandom(0.003,1);
         int s;
         for (s = 2; s < DimPop; s++) {
             if (SelectionWheel[s] < r && r <= SelectionWheel[s + 1]) {
                 break;
             }
         }
          population[i]=population[s];
    }
}

void TournamentSelection(vector<Cromozom> &population,int k)
{

   
   int d=population[0].route_dimension;
   float min=1000;
   
    for(int i=2;i<DimPop;i++)
     {
        Cromozom best=population[i];
        int l = k;  // asta e la misto
       
        int v[d]={0}; 
       
         for(int j=0;j<l;j++)
         {
            int a=generateRandom(0,d-1);
            //while(v[a] != 0)   // de gandit!
             //   {
              //      int a=generateRandom(0,d-1);
               // }
            if(population[a].getFitness()<min)
                  {
                    min=population[a].getFitness();
                    best=population[a];
                  }
            v[a]=1;       
         }
         population[i]=best;
     }
}


Round_RobinTournamentSelection(vector<Cromozom> &population,int K) /*RTS*/
{
  //!!!Mai jos eu am impartit in cele 2 grupuri (prima jum ; a 2 a jumatate) dar as dori sa alegem elemntele random in cele 2 parti ;

    float TotalScore=0;
    float TotalScore2=0;
    float Distrib[300];
    float Distrib2[300];
    float ProbSurv[300];
           
    for(int i=1;i<DimPop/2;i++)
         TotalScore += population[i].getFitness();
 
    for(int i=DimPop/2+1;i<DimPop;i++)
         TotalScore2 += population[i - 1].getFitness(); //acel omega 0

    for(int i=1;i<DimPop/2;i++)
         Distrib[i] = population[i - 1].fitness/TotalScore; //acel omega 0 ; probabilitatea de alegere in grupul 1
    
    for(int i=DimPop/2+1;i<DimPop;i++)
         Distrib2[i-DimPop] = population[i - 1].fitness/TotalScore2; //acel omega 0 ; probabilitatea de alegere in grupul 2
        
         float Sum_j=0;
        
        for(int i=DimPop/2;i<DimPop;i++)
         Sum_j+=K-(i-DimPop);

   //Sum(Distrib[i])=1; cam asta e regula
    
    for(int i=1;i<DimPop/2;i++)
         ProbSurv[i]=2*(i-1)/(K*(K-1))*Distrib[i] + 4*(1-Distrib[i])/((K*K)*(K-1))*Sum_j;//vezi la acel i ca ia val dubioase
         //Si da in final nu avem 2 valori cu aceeasi probabilitate/
        
     for(int i=DimPop/2+1;i<DimPop;i++)
        ProbSurv[i]=2*(i-1)/(K*(K-1))*Distrib[i] + 4*(1-Distrib[i])/((K*K)*(K-1))*Sum_j;//vezi la acel i ca ia val dubioase
         //Si da in final nu avem 2 valori cu aceeasi probabilitate/

   int d=population[0].route_dimension;
   float min=1000;
   
    for(int i=1;i<DimPop;i++)   //Nu conteaza grupul ; se alege random un elemt ;
     {
        Cromozom best=population[i];
        int l = 48;  // asta e la misto pt pool size
       
        int v[d]={0}; 
       
         for(int j=0;j<l;j++)
        {
            int a=generateRandom(0,d-1);
         /*   float T=100;
            while(v[a]!=0 && T>1)
            {  
              if( generateFloatRandom(0,1-0.00001)< exp(abs(ProbSurv[i]-min)/ T ))    
              {
                {
                    v[a]=0;
                    T=T*0.1;
                }
              }
              else
               {
                 int a=generateRandom(0,d-1);
               }
            }
         */

            //while(v[a] != 0)   // de gandit!
             //   {
              //      int a=generateRandom(0,d-1); //as vrea sa iau mereu in pool elemente diferite cat de cat(functie de tem facem in cazul in care se nim acelasi elemnt)
               // }
            if(ProbSurv[i]<min)
                  {
                    min=ProbSurv[i];
                    best=population[a];
                  }
            
            v[a]=1;       
        }


         //Pressure : Cumva as vrea cam maxim 2 sa fie la fel in populatie (iar la o dublura sau dau asa o sansa de reinsertie dupa presiune)
     /*    float Temp=1;
        
        for(int j=0;j<i;j++)
           if(population[j] == best)
              Temp=Temp*0.5;

          float a = generateFloatRandom(0,1-0.00001);

           if(a<Temp)
            population[i]=best;
            else
           i--; 
        */
         population[i]=best;   
     }

  //!!FINAL COMPLEXITY O(N^2) fata de restul cu max NlogN sau N in maj cazurilor :))
}

/*

Basic types of tournaments.
Tournament formats :

1. Single-elimination tournament.        Done it ;
2. Double-elimination tournament.        Done it ;
3. Round Robin.                          Done it;
4. Snake (a partial Round Robin)         In progress;


/*

Here we try to use a HyperCube   Ce a mai ramas;

void SnakeTournament(vector<Cromozom> &population)
{

}
*/



void ClusteringSelection(vector<Cromozom> &population) //merge cu turneu selection
{
    sort(population.begin(), population.end(), orderByFitness);
    
    int NumCluster=0;
    float values[201];
    int Hash[201];           //Multime in care fiecare element al pop are valoare clusterului din care face parte;
    float Sum=0;

    for(int i=1;i<DimPop;i++)
    {
        values[i - 1]=population[i].getFitness();
        Sum=Sum+values[i];
    }

  //Aici calculam constantele pentru formarea clusterelor;  
    float mean=Sum/199;
    float sd=standardDeviation(values,199);
    float alfa=sqrt(sd);

   //Facem impartirea in Clustere : (salvez si nr-ul lor)
      
      for(int i=1;i<DimPop-1;i++)
       {
             NumCluster++;
             int ClusterDim[200];// Sper ca e bine acum ; am corectat if-ul ;
             for(int j=i+1;j<DimPop;j++)
                if(abs(population[i].getFitness() - population[j].getFitness())<=alfa && ( (Hash[i]==0 && Hash[j]==0) or ( ClusterDim[i] > 1 && Hash[i] == 1 && Hash[j] == 0 ) or (ClusterDim[j] > 1 && Hash[i] == 0 && Hash[j] == 1 ) ) )
                   {
                       Hash[i]=NumCluster;
                       Hash[j]=NumCluster;
                       ClusterDim[i]++;
                       ClusterDim[j]++;
                   }
       }                                                  

   //Evaluam Constantele pentru Clustere : SD_Cluster ; CV_CLuster ; MEAN_CLuster ; Cluster_dimension ; difference D = sqrt(SD_Cluster)
      //Start Cluster Reodronation.   (o voi face prin reatribuirea permutarilor (vad 1 pe ce loc e de fapt 2 .... n pe ce loc este )) prin vect de atribuiri v[1]=k(in reord 1 e pe poz k)

 float MatCluster[200][DimPop];
 int ClusterDim[200];
 float Mean_Cluster[200];
 float SD_Cluster[200];
 float D_Cluster[200];
 float CV_Cluster[200];

 for(int i=0;i<NumCluster;i++)
 {  
    int k=0;
    ClusterDim[i]=0;
    int S=0;
    float valuess[200];
    for(int j=0;j<DimPop;j++)
    {
        if(Hash[j]==i)
       {
        MatCluster[i][k]=values[j];
        k++;
        ClusterDim[i]++;
        S+=values[j];
        valuess[j]=values[j];
       }
         else
      {
        MatCluster[i][k]=0;
        k++;
      }
    }

    Mean_Cluster[i]=S/ClusterDim[i];
    SD_Cluster[i]=standardDeviation(valuess,ClusterDim[i]);
    D_Cluster[i]=sqrt(SD_Cluster[i]);
    CV_Cluster[i]=(SD_Cluster[i]*ClusterDim[i])/(D_Cluster[i]+1); //Acum putem in sfarsit sa le ordonam dupa CV_Cluster.  //D_Cluster[i]*MatCluster[i][]+1
 }
 
int ReordCl[200];
int RevReord[200];
 for(int i=0;i<NumCluster;i++)
{
    int k=1;
  for(int j=0;j<NumCluster;j++)
    if( i!=j && CV_Cluster[i]>CV_Cluster[j])   //< nu mai stiu care trb;
            k++;
  ReordCl[i+1]=k; 
  RevReord[k]=i+1;               //Cu astea reordonez populatia dupa clusterul maxim din care face parte    
}       
 
//Refac Hash-ul

vector<Cromozom> CopyPop; //Aici copia pop pe care voi face operatiile

for(int i=1;i<DimPop;i++)
{
    for(int j=1;j<=NumCluster;j++)
    CopyPop.emplace_back(population[RevReord[j]]); //Sper ca e oke
}

  vector<Cromozom> NewPop(300);
 
 int k=1; //De la 1 pt ca avem elitism pe prima valoare
   //Parintii vor fi din acelasi cluster - cei mai buni 
   //Nu sunt sigur dar vreau ca sa i pun pe perechi crescator (ma gandesc ca asa vor face crossOver)
  
    for(int i=1;i<=NumCluster;i++)            ///
    {
        int sw=0;
        for(int j=1;j<DimPop;j++)   //Am ales primele 2 elemente dar inainte am ordonat dupa constanta.
        {
          if(Hash[j]==i && sw==0)
          {
            sw++;
            NewPop[k]=CopyPop[j];   //Aici cum sa af=daug pe pozitia k CopyPop[j]
            k++;
          }
         else
         {
           if(Hash[j]==i && sw==1)
           {
            sw++;
            NewPop[k]=CopyPop[j];     //La fel
            k++;
           }
           else
           {
              if(sw==2)
                sw=0;     //Resetez. Am facut asa ca poate e nev sa le controlez cumva.
           }
         } 
        }
    }
    
   //Aici updatez populatia mea;
    for(int i=1;i<DimPop;i++)             ///Aceeasi problema
      population[i]=NewPop[i];
     

  //   for(int i=1;i<DimPop;i++)
       //    population[i].emplace_back(NewPop[i]);
  //  population.emplace_back(NewPop);
}


void Evaluate(vector<Cromozom> &population)
{
        for(auto iter = population.begin(); iter != population.end(); iter++)
            iter->update_fitness();
            
        sort(population.begin(), population.end(), orderByFitness);
}

void PopulationGeneration(vector<Cromozom> &PopulationA)
{   
    Cromozom init_crom = Read("input/150.tsp");
   // init_crom.apply_SA();
    PopulationA.emplace_back(init_crom);
    for(int i=0;i<DimPop;i++)
    {   
        init_crom.random_route();
      //  init_crom.apply_SA();
        PopulationA.emplace_back(init_crom);
    }   
}

void print(const vector<Cromozom> &pop, int size)
{
    for(auto iter = pop.begin(); iter != pop.end(); ++iter)
        for(int j = 0; j < iter->route_dimension; j++)
            cout<<"INDEX: " << j <<" X:" << iter->route[j].x << "Y: " << iter->route[j].y << endl;
}
