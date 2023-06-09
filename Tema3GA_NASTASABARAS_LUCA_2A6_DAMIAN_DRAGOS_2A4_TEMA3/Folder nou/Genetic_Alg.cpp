#include<iostream>
#include "Operators.cpp"

#define MAX_GENERATION 2000


void SimulatedEvolution()   
{
    ofstream out;
    char* filename = "Rezultate1.txt";
    out.open(filename);
    vector<Cromozom> population;
    Cromozom curr_best {};
    PopulationGeneration(population);
    Evaluate(population);
    cout << "BEST INIT COST: " << population[0].getDistance() << endl;
    curr_best = population[0];
    for(int gen = 1; gen <= MAX_GENERATION; gen++)
    {
       // Selectie:
       
       //SelectWheelOfFortune(population);
       // TournamentSelection(population,18);
       ClusteringSelection(population);
       // Round_RobinTournamentSelection(population,200);
       // cout<<"iese"<<'\n';
     
       //CrossOver:

        IGXCrossOver(population);
       // CXCrossOver(population);
       // ElitistPopCross1(population);
       //OXCrossover(population);
       
       //Mutation:
       
       // Mutate(population);
       // Inverse(population);
       // RGIBNM(population);
       // IRGIBNM(population);
       SBM(population);//Have the best from those 4 mutation
        //     HyperMutate(population);
       
       //Evaluate:
        Evaluate(population);
        out<<population[0].getDistance()<<'\n';  //-1900
    }

    cout << "BEST COST AG: " << population[0].getDistance() << endl;
    cout << "Ordinea de vizitare este: ";
    for (auto const &city : population[0].route)
        cout << city.index << " ";  
    cout <<endl;
}