#include<iostream>
#include"Cromozom.cpp"
#define SimulatedAnnealingIterrations 10000

void SA()
{
    Cromozom a;
    Cromozom neigbour;
    float finalCost = 0;
   
    a=Read("input/52.tsp");
    a.random_route();
   ofstream out;
    char* filename = "Rezultate1.txt";
    out.open(filename);
    float Cost1=0;  
    float Cost2=0;

   for(int i=1;i<=SimulatedAnnealingIterrations;i++)
   {
      Cost1=a.getFitness();     
      float t=100;
     
      while(t>0.000001)
      {
         neigbour=a; 
       //  neigbour.swapp_hemming();                     
         neigbour.swapp_reverse();
         Cost2=neigbour.getFitness(); 
         if(Cost1>Cost2)   //<
         {
            a=neigbour;
            Cost1=Cost2;
         }
         else
         {
            float abso = abs(Cost2 - Cost1);
            float exponent = exp( -abso / t );
            float random = generateFloatRandom(0, 1 - 0.000001);
            if(random < exponent)
            {
               a = neigbour;
               Cost1=Cost2;
               t=t*0.9995;
            }
         }
      }
      Cost1 = max(Cost1, Cost2);//min

      finalCost=a.getDistance();
      out<<finalCost<<'\n';
      
      cout << "ITER: " << i << endl;
   }
   
  // cout << "FINAL " << finalCost << endl;
 //  for(int i=0;i<a.route_dimension;i++)
 //  cout<<a.route[i]<<' ';
}