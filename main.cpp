#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <assert.h>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>
#include <sstream>


/*
SORUNLAR
1)
*/

using namespace std;

//float min = 1000;
float minimum = 1000;
float maximum = 0;



const int pop_Size = 10;
const float KE_Loss_Rate = 0.2;
const float MoleColl = 0.4;
int buffer = 0;
int Initial_KE = 1000;
const int alpha = 5;
const int beta = 15000;

const int num_node = 500;
const int num_molecule = 10;
const int num_parents = 5;
const int range = 800;
const int num_facility = 5;
const int iter = 5;





const float mutation_rate = 0.3;
const float replacement_prob = 0.3;
const int num_iter = 20;
const float prob_repairing = 0.5;
string problem_demand = "sjc402_demand.txt";
string problem_distance = "sjc402_distance.txt";
string problem_result = "sjc402-1p-max_imp.txt";
string improving="a";
string xover_operator="b";

struct node {
	double x_cord;
	double y_cord;
	int demand;
	int satisfied_demand;
};

node Array_node[num_node];
float distance_node[num_node][num_node];

struct molecule {
	int soln[num_node];
	int satisfied[num_node];
	int PE;///Potential Energy = fitness value
	int KE;///Kinetic Energy
	int NumHit;
	int MinStruct[num_node];
	int MinPE;
	int MinHit;
	float prob;
	int index;
};

struct population {
	molecule chros[num_molecule];
};

struct parent {
	///molecule parents[num_parents];
	vector<molecule> parents[num_parents];
	///molecule<vector> parents[num_parents];
};

vector <molecule> molecules_vec;

molecule m1;
molecule m2;
molecule ow_m1;
molecule syn_m;
molecule im_m1;
molecule im_m2;

parent mating_pool;

population initial;

void init_pop(void)  ///Generating Initial Population
{
	for (int i = 0; i < num_molecule; i++)
    {
		int counter = 0;
		while (counter <= num_facility - 1)
		{
			int rand_num = rand() % num_node;
			if (initial.chros[i].soln[rand_num] == 0)
            {
				initial.chros[i].soln[rand_num] = 1;
				counter += 1;
			}
		}
		initial.chros[i].index = i;
        initial.chros[i].KE = Initial_KE;
        initial.chros[i].NumHit = 0;
        initial.chros[i].MinHit = 0;
	}
}

int total_demand[num_node];

void demand_satisfied_node() { ///Calculates each nodes' satisfaction portion in case facility opened

	for (int i = 0; i < num_node; i++)
	{
		for (int j = 0; j < num_node; j++)
		{
			if (distance_node[i][j] < range) {
				total_demand[i] += Array_node[j].demand;
			}
		}
		Array_node[i].satisfied_demand = total_demand[i] + Array_node[i].demand;
	}
}

void ordering() { ///Ordering nodes according to their demand_satisfied_node()
	int i, j;
	node temp;

	for (i = 0; i < num_node; i++)
	{
		for (j = i + 1; j < num_node; j++)
		{
			if (Array_node[i].satisfied_demand < Array_node[j].satisfied_demand)
			{
				temp = Array_node[i];
				Array_node[i] = Array_node[j];
				Array_node[j] = temp;
				i = 0;
			}
		}
	}
}

void satisfied_vec_func(molecule &x)///satisfied vector calculation
{
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}

	for (int i = 0; i < num_node; i++)
    {
		for (int j = 0; j < num_node; j++)
        {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range))
            {
				x.satisfied[j] = 1;
			}
			else
            {
				if(x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
}

void satisfied_vec_func_offspring(molecule &x) { //satisfied vector calculation for offsprings
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}
	for (int i = 0; i < num_node; i++) {

		for (int j = 0; j < num_node; j++) {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range)) {
				x.satisfied[j] = 1;
			}
			else {
				if (x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
}

int fitness_func(molecule x) { //Fitness function calculation
	x.PE = 0;
	int sum = 0;
	for (int i = 0; i < num_node; i++) {
		sum = x.satisfied[i] * Array_node[i].demand + sum;
		/*
		if(x.satisfied[i]==1)
        {
            cout << Array_node[i].demand << " ";///kon kon demand value add kore fitness value passi seta
        }
        */
	}
	///cout << "Sum=";
	///cout << sum << endl;
	x.PE = sum;
	return x.PE;
}

void prob_func(population &gen) /// Calculating probabilities of chromosomes ***Note: Min Fitness Prob=0***
{
	for (int i = 0; i < num_molecule; i++)
    {
		if (gen.chros[i].PE < minimum)
        {
			minimum = gen.chros[i].PE;
		}
		else if (gen.chros[i].PE > maximum)
		{
			maximum = gen.chros[i].PE;
		}
	}
	///cout << maximum << endl;
	///cout << minimum << endl;
	for (int i = 0; i < num_molecule; i++) {
		gen.chros[i].prob = (gen.chros[i].PE - minimum) / (maximum - minimum);
	}
	float sum = 0;
	for (int i = 0; i < num_molecule; i++)
	{
		sum += gen.chros[i].prob;
	}
	///cout << sum;
	gen.chros[0].prob = gen.chros[0].prob / sum;
	for (int i = 1; i < num_molecule; i++)
	{
		gen.chros[i].prob = gen.chros[i - 1].prob + (gen.chros[i].prob / sum); //Finding Cumulative Probabilities
	}
}

float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

/*
void parent_selection(population gen) { ///Selecting # of distinct parents out of population

	int m = -1;
	int iteration = 0;

	for (int h = 0; h < num_parents + iteration; h++)
	{
		m++;
		float rnd = RandomFloat(0, 1);
		float diff = 1;
		for (int j = 0; j < num_molecule; j++)
		{
			if (((gen.chros[j].prob - rnd) > 0) && ((gen.chros[j].prob - rnd) < diff)) {
				diff = abs(gen.chros[j].prob - rnd);
                mating_pool.parents[m] = gen.chros[j];
                ///mating_pool.parents.push_back(gen.chros[j]);
				mating_pool.parents[m].index = gen.chros[j].index;
				///mating_pool.parents[m].index = gen.chros[j].index;


				for (int k = 0; k < m; k++)
				{
					if (mating_pool.parents[k].index == mating_pool.parents[m].index)
					{
						m = m - 1;
						iteration++;
					}
				}
			}
		}
	}
}

*/

molecule crossover(molecule parent1, molecule parent2) //Crossover Operator(1-point)
{
	xover_operator = "1-Point Crossover";
	int k = 188; //Cut Point
	molecule offspring;
	for (int i = 0; i < k; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k; i < num_node; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	return offspring;
}

molecule crossover_2p(molecule parent1, molecule parent2) //Crossover Operator(2-point)
{
	xover_operator = "2-Point Crossover";
	int k1 = 180; //Cut Point
	int k2 = 246;
	molecule offspring;
	for (int i = 0; i < k1; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	for (int i = k1; i < k2; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k2; i < num_node; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	return offspring;
}

molecule tri_crossover(molecule parent1, molecule parent2, molecule parent3) //3-Parent Crossover
{
	int k1 = 182;
	int k2 = 265;
	molecule offspring;
	for (int i = 0; i < k1; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	for (int i =k1; i <k2; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k2; i < num_node; i++)
	{
		offspring.soln[i] = parent3.soln[i];
	}
	return offspring;
}

molecule satisfied_new(molecule x) { //satisfied vector calculation
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}
	for (int i = 0; i < num_node; i++) {

		for (int j = 0; j < num_node; j++) {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range)) {
				x.satisfied[j] = 1;
			}
			else {
				if (x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
	return x;
}

molecule max_imp_repairing(molecule x) {
	improving = "Maximum Improving";
	int sum = 0;
	int max_demand = 0;
	int min_demand = 10000;
	int index = 0;
	x = satisfied_new(x);
	x.PE = fitness_func(x);
	int current_fitness = x.PE;
	molecule y;
	y = x;
	int deviation[num_node];

	for (int i = 0; i < num_node; i++)
	{
		sum += x.soln[i];
	}
	if (sum < num_facility) {
		while (sum < num_facility)
		{
			for (int i = 0; i < num_node; i++)
			{
				if (y.soln[i] == 0)
				{
					y.soln[i] = 1;
					y=satisfied_new(y);
					y.PE=fitness_func(y);
					deviation[i] = y.PE - current_fitness;
				}
				y = x;
			}
			for (int i = 0; i < num_node; i++)
			{
				if (deviation[i] > max_demand)
				{
					max_demand = deviation[i];
					index = i;
				}
			}

			x.soln[index] = 1;
			sum++;
		}
	}

	if (sum > num_facility) {
		while (sum > num_facility)
		{
			for (int i = 0; i < num_node; i++)
			{
				if (y.soln[i] == 1)
				{
					y.soln[i] = 0;
					y=satisfied_new(y);
					y.PE = fitness_func(y);
					deviation[i] = current_fitness-y.PE;
				}
				y = x;
			}
			for (int i = 0; i < num_node; i++)
			{
				if ((deviation[i] <min_demand)&&(deviation[i]>0))
				{
					min_demand = deviation[i];
					index = i;
				}
			}

			x.soln[index] = 0;
			sum--;
		}
	}
	return x;
}

molecule repairing(molecule x) { //Randomly Repairing
	improving = "Random Repairing";
	int sum = 0;
	for (int i = 0; i < num_node; i++)
	{
		sum += x.soln[i];
	}

	if (sum < num_facility) {
		while (sum <= num_facility - 1)
		{
			int rand_num = rand() % num_node;
			if (x.soln[rand_num] == 0) {
				x.soln[rand_num] = 1;
				sum += 1;
			}
		}
	}
	if (sum > num_facility) {
		while (sum >= num_facility + 1) {
			int rand_num = rand() % num_node;
			if (x.soln[rand_num] == 1) {
				x.soln[rand_num] = 0;
				sum -= 1;
			}
		}
	}
	return x;
}

molecule alt_repairing(molecule x) { //Repairing according to demand_satisfied_node
	improving = "Greedy Repairing";
	int sum = 0;
	float rand_num = 0;
	for (int i = 0; i < num_node; i++)
	{

		sum += x.soln[i];
	}

	int maks = 0;

	while (sum < num_facility) {
		for (int i = 0; i < num_node; i++)
		{
			rand_num = RandomFloat(1, 0);
			if ((x.soln[i] == 0) && (rand_num < prob_repairing))
			{
				x.soln[i] = 1;
				sum++;
			}
		}
	}

	while (sum > num_facility)
	{
		for (int i = num_node; i > -1; i--)
		{
			rand_num = RandomFloat(1, 0);

			if ((x.soln[i] == 1) && (rand_num < prob_repairing))
			{
				x.soln[i] = 0;
				sum--;
			}
		}
	}

	return x;
}

molecule mutation(molecule x) { ///Mutation operator (swap two elements of array)
	float rand_num = RandomFloat(1, 0);
	int rnd = rand() % num_node;
	if (rand_num < mutation_rate) {
		for (int i = 0; i < num_node; i++)
		{
			swap(x.soln[rnd], x.soln[rand() % num_node]);
		}
	}
	return x;
}

void inter_molecular_ineffective_collision(molecule m1, molecule m2, int current_index, int next_index)
{
    int copy_index = num_facility-1;
    int break_index;
    int q = 0;
    int c = 0;
    int r;
    int p;
    int E_inter;
    float rf;

    ///creation of first molecule from m1 molecule
    for(int i = 0; i < num_node; i++)
    {
        if(m1.soln[i]==1)
        {
            c++;
        }
        if(c <= copy_index)
        {
            im_m1.soln[i] = m1.soln[i];
        }
        if(c > copy_index)
        {
            im_m1.soln[i] = 0;
        }
        if(c==copy_index && q==0)
        {
            break_index = i;
            q++;
        }
    }

    r = num_node - break_index;
    p = rand() % r + break_index;
    im_m1.soln[p] = 1;
    ///first molecule creation finished


    ///creation of second molecule from m2 molecule
    c = 0;
    q = 0;
    for(int i = 0; i < num_node; i++)
    {
        if(m2.soln[i]==1)
        {
            c++;
        }
        if(c <= copy_index)
        {
            im_m2.soln[i] = m2.soln[i];
        }
        if(c > copy_index)
        {
            im_m2.soln[i] = 0;
        }
        if(c==copy_index && q==0)
        {
            break_index = i;
            q++;
        }
    }

    r = num_node - break_index;
    p = rand() % r + break_index;
    im_m2.soln[p] = 1;
    ///creation of second molecule finished

    /*
    cout << "m1 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m1.soln[i];
    }
    cout << endl;

    cout << "im_m1 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << im_m1.soln[i];
    }
    cout << endl;


    cout << "m2 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m2.soln[i];
    }
    cout << endl;

    cout << "im_m2 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << im_m2.soln[i];
    }
    cout << endl;
    */

    satisfied_vec_func(im_m1);
    im_m1.PE = fitness_func(im_m1);
    satisfied_vec_func(im_m2);
    im_m2.PE = fitness_func(im_m2);
    /*
    cout << "m1.PE = "<< m1.PE << endl;
    cout << "im_m1.PE = "<< im_m1.PE << endl;
    cout << "m2.PE = "<< m2.PE << endl;
    cout << "im_m2.PE = " <<im_m2.PE << endl;
    */
    molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;//current_index++;
    molecules_vec[next_index].NumHit = molecules_vec[next_index].NumHit + 1;

    E_inter = (im_m1.PE + im_m2.PE) - (m1.PE + m2.PE);

    if(E_inter >= 0)
    {
        rf = RandomFloat(0,1);
        im_m1.KE = E_inter * rf;
        im_m2.KE = E_inter * (1 - rf);
        molecules_vec[current_index].PE = im_m1.PE;
        molecules_vec[next_index].PE = im_m2.PE;
        molecules_vec[current_index].KE = im_m1.KE;
        molecules_vec[next_index].KE = im_m2.KE;

        if(im_m1.PE > molecules_vec[current_index].MinPE)
        {
            for(int i = 0; i < num_node; i++)
            {
                molecules_vec[current_index].MinStruct[i] = im_m1.soln[i];
            }
            molecules_vec[current_index].MinPE = im_m1.PE;
            molecules_vec[current_index].MinHit = molecules_vec[current_index].NumHit;
        }


        if(im_m2.PE > molecules_vec[next_index].MinPE)
        {
            for(int i = 0; i < num_node; i++)
            {
                molecules_vec[next_index].MinStruct[i] = im_m2.soln[i];
            }
            molecules_vec[next_index].MinPE = im_m2.PE;
            molecules_vec[next_index].MinHit = molecules_vec[next_index].NumHit;
        }
    }


}

void synthesis(molecule m1, molecule m2, int current_index, int next_index)
{
    int c = 0;
    int desire_index = num_facility - 1;
    int target_index;

    for(int i = 0; i < num_node; i++)
    {
        if(m1.soln[i]==1)
        {
            c++;
        }
        if(c <= desire_index)
        {
            syn_m.soln[i] = m1.soln[i];
        }

        if(c == desire_index)
        {
            target_index = i;
            break;
        }
    }


    for(int i = target_index+1; i < num_node; i++)
    {
        if(m2.soln[i]==1)
        {
            c++;
        }
        if(c <= num_facility)
        {
            syn_m.soln[i] = m2.soln[i];
        }
        if(c > num_facility)
        {
            syn_m.soln[i] = 0;
        }
    }

    /*
    cout << "m1 molecule" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m1.soln[i];
    }
    cout << endl;

    cout << "m2 molecule" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m2.soln[i];
    }
    cout << endl;

    cout << "syn_m molecule" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << syn_m.soln[i];
    }
    cout << endl;
    */
    satisfied_vec_func(syn_m);
    syn_m.PE = fitness_func(syn_m);
    /*
    cout << m1.PE << endl;
    cout << m2.PE << endl;
    cout << syn_m.PE << endl;
    */
    ///if((m1.PE + m2.PE + m1.KE + m2.KE) >= (syn_m.PE))
    if((syn_m.PE >= m1.PE) || (syn_m.PE >= m2.PE))
    {
        syn_m.KE = (m1.PE + m2.PE) - syn_m.PE;
        for(int i = 0; i < num_node; i++)
        {
            syn_m.MinStruct[i] = syn_m.soln[i];
        }
        syn_m.MinPE = syn_m.PE;

        if(m1.PE < syn_m.PE)
        {
            molecules_vec.erase(molecules_vec.begin() + current_index);
        }
        if(m2.PE < syn_m.PE)
        {//current_index++;
            molecules_vec.erase(molecules_vec.begin() + next_index);
        }

        molecules_vec.push_back(syn_m);
    }
    else
    {
        molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;
        //current_index++;
        molecules_vec[next_index].NumHit = molecules_vec[next_index].NumHit + 1;
    }

}

void on_wall(molecule m, int current_index)
{
    int one_count = num_facility - 1;
    int desire_index;
    int r;
    int p;
    int c = 0;
    float a;


    ///creation of new solution from m molecule
    for(int i = 0; i < num_node; i++)
    {
        if(m.soln[i]==1)
        {
            c++;
        }
        if(c <= one_count)
        {
            ow_m1.soln[i] = m.soln[i];
        }
        if(c > one_count)
        {
            ow_m1.soln[i] = 0;
        }
        if(c == one_count)
        {
            desire_index = i;
        }

    }
    r = num_node - desire_index;
    p = rand() % r + desire_index;
    ow_m1.soln[p] = 1;
    ///molecule creation done

    satisfied_vec_func(ow_m1);
    ow_m1.PE = fitness_func(ow_m1);

    /*
    cout << "This is m:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m.soln[i];
    }
    cout << "\t";
    cout << m.PE <<endl;

    cout << "This is ow_m1:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << ow_m1.soln[i];
    }
    cout << "\t";
    cout << ow_m1.PE;
    */
    molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;

    ///if((m.PE + m.KE) >= ow_m1.PE)
    if(ow_m1.PE >= m.PE)
    {
        a = RandomFloat(KE_Loss_Rate,1);
        ///ow_m1.KE = (m.PE - ow_m1.PE + m.KE)*a;
        ow_m1.KE = abs(ow_m1.PE - (m.PE + m.KE))*a;
        ///buffer = buffer + ((m.PE - ow_m1.PE + m.KE)*(1 - a));
        buffer = buffer + (abs(ow_m1.PE - (m.PE + m.KE)))*(1 - a);

        for(int i = 0; i < num_node; i++)
        {
            molecules_vec[current_index].soln[i] = ow_m1.soln[i];
        }
        molecules_vec[current_index].PE = ow_m1.PE;
        molecules_vec[current_index].KE = ow_m1.KE;

        if(molecules_vec[current_index].PE > molecules_vec[current_index].MinPE)
        {
            for(int i = 0; i < num_node; i++)
            {
                molecules_vec[current_index].MinStruct[i] = molecules_vec[current_index].soln[i];
            }
            molecules_vec[current_index].MinPE = molecules_vec[current_index].PE;
            molecules_vec[current_index].MinHit = molecules_vec[current_index].NumHit;
        }
    }
}


void decomposition(molecule m, int current_index)
{
    int remain;
    int q;
    int r;
    int p;
    int j;
    int c = 0;///1 counter in molecule
    int E_dec;
    float rf1;
    float rf2;
    float rf3;
    float rf4;
    //molecule m1;
    //molecule m2;

    int k = num_node/2 - 2;
    int k1 = k;


    ///creation of m1 from m
    for (int i = 0; i < num_node; i++)
    {
        if(i < k)
        {
            if(m.soln[i] == 1)
            {
                c++;
            }
            m1.soln[i] = m.soln[i];
        }
        else{
            m1.soln[i] = 0;
        }
    }

    remain = num_facility - c;
    for(int i = 0; i < remain; i++)
    {
        r = num_node - k1;
        p = rand() % r + k1;
        m1.soln[p] = 1;
        k1++;
    }


    ///creation of m2 from m
    j = 0;
    c = 0;
    for(int i = k; i < num_node; i++)
    {
        if(m.soln[i] == 1)
        {
            c++;
        }
        m2.soln[j] = m.soln[i];
        j++;
    }

    for(int i = j; i < num_node; i++)
    {
        m2.soln[i] = 0;
    }

    remain = num_facility - c;
    for(int i = 0; i < remain; i++)
    {
        r = num_node - k;
        p = rand() % r + k;
        m2.soln[p] = 1;
        k++;
    }

    /*
    cout << "This is m:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m.soln[i];
    }

    cout << "This is m1:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m1.soln[i];
    }

    cout << "This is m2:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m2.soln[i];
    }
    */
    satisfied_vec_func(m1);
    m1.PE = fitness_func(m1);
    satisfied_vec_func(m2);
    m2.PE = fitness_func(m2);

    /*
    cout << m.PE << endl;
    cout << m1.PE << endl;
    cout << m2.PE << endl;
    */
    ///if((m.PE + m.KE) >= (m1.PE + m2.PE))
    if((m1.PE >= m.PE) || (m2.PE >= m.PE))
    {
        E_dec = (m1.PE + m2.PE) - (m.PE + m.KE);
        ///cout << "E_dec = " << E_dec;
        rf1 = RandomFloat(0,1);
        ///cout << " rf1 = " << rf1;
        m1.KE = E_dec * rf1;
        m2.KE = E_dec * (1-rf1);

        ///cout << " m1.KE = " << m1.KE;
        ///cout << " m2.KE = " << m2.KE << endl;



        for(int i = 0; i < num_node; i++)
        {
            m1.MinStruct[i] = m1.soln[i];
            m2.MinStruct[i] = m2.soln[i];
        }
        m1.MinPE = m1.PE;
        m2.MinPE = m2.PE;

        molecules_vec.erase(molecules_vec.begin() + current_index);

        molecules_vec.push_back(m1);
        molecules_vec.push_back(m2);

    }
    else
    {
        rf2 = RandomFloat(0,1);
        rf3 = RandomFloat(0,1);

        ///E_dec = (m.PE + m.KE + ((rf2 * rf3) * buffer)) - (m1.PE + m2.PE);
        E_dec = (m1.PE + m2.PE + ((rf2 * rf3) * buffer)) - (m.PE + m.KE);

        if(E_dec >= 0)
        {
            buffer = buffer * (1 - (rf2 * rf3));
            rf4 = RandomFloat(0,1);
            m1.KE = E_dec * rf4;
            m2.KE = E_dec * (1 - rf4);

            for(int i = 0; i < num_node; i++)
            {
                m1.MinStruct[i] = m1.soln[i];
                m2.MinStruct[i] = m2.soln[i];
            }

            m1.MinPE = m1.PE;
            m2.MinPE = m2.PE;

            if((m1.PE > m.PE) || (m2.PE > m.PE))
            {
                molecules_vec.erase(molecules_vec.begin() + current_index);
            }
            if(m1.PE > m.PE)
            {
                molecules_vec.push_back(m1);
            }
            if(m2.PE > m.PE)
            {
                molecules_vec.push_back(m2);
            }
        }
        else
        {
            molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;
        }
    }
         /*
         for (int i = 0; i < molecules_vec.size(); i++)
         {
             for(int j = 0; j < num_node; j++)
             {
                 cout << molecules_vec[i].satisfied[j];
             }
             cout << "\t";
             cout << molecules_vec[i].PE << endl;
         }
         */
}

int main()
{
	int bekle;
	int mak_out = 0;
	ofstream result("Result.txt");
	ofstream test1("single_facility_all_result_sjc500_800.txt");
	double timer=0;
	double timer_out = 0;
	int sum_time = 0;
	int sum_best = 0;
	int ran;
	int ran1;
	int ran2;
	int first;
	int second;
    int best_single = 0;
    int best_single_index;
    int best_second_index;
    int best_third_index;
    int best_fourth_index;
    int ind;
    int value = 0;
    int ind1;
    int value1 = 0;
    int ind2;
    int value2 = 0;
    int data[num_node];

		int start_s = clock();
		// Start timer

		srand((int)time(NULL));

		ofstream out("sjc324-1p-max_imp.txt");
		ifstream get_demand; //Get demand values from problem set

		get_demand.open("sjc500_demand.txt");
		for (int i = 0; i < num_node; i++)
		{
			get_demand >> Array_node[i].demand;
		}

		get_demand.close();

		ifstream get_distance; //Get distance values from problem set

		get_distance.open("sjc500_distance.txt");

		for (int i = 0; i < num_node; i++)
		{
			for (int j = 0; j < num_node; j++)
			{
				get_distance >> distance_node[i][j];
			}
		}

		get_distance.close();

		demand_satisfied_node(); //Calculate each node's portion of satisfaction

		//ordering(); //Order nodes according to portion of satisfaction

        /*
		molecule x;
		for(int j = 0; j < num_node; j++)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==j)
                {
                   x.soln[i] = 1;
                }
                else{
                   x.soln[i] = 0;
                }
            }

            for(int i = 0; i < num_node; i++)
            {
               cout << x.soln[i];
            }
            cout << "\n";
            cout << j << " ";
            test1 << j << " ";
            satisfied_vec_func(x);
            x.PE = fitness_func(x);
            if(x.PE > best_single)
            {
                best_single = x.PE;
                best_single_index = j;
            }
            cout << x.PE << endl;
            test1 << x.PE << endl;
            for(int p = 0; p < num_node; p++)
            {
                cout << x.satisfied[p];
            }
            cout << "\n\n";

        }
        test1 << "\n\n";
        test1 << best_single_index << " number index is 1 among 500 indices when get the best solution" << endl;
        test1 << "Result : " << best_single;
        cout << best_single_index << " ";
        cout << "Result : " << best_single <<endl;
        */





         ///total second
         /*
         molecule y;
         int u = 0;
         while(u < 500)
         {
            for(int j = 0; j < num_node; j++)
            {
                for(int i = 0; i < num_node; i++)
                {
                   if(i==u || i==j)
                   {
                      y.soln[i] = 1;
                   }
                   else
                   {
                      y.soln[i] = 0;
                   }
                }
                ///for(int i = 0; i < num_node; i++)
                ///{
                ///   cout << y.soln[i];
                ///}
                ///cout << "\n";
                satisfied_vec_func(y);
                y.PE = fitness_func(y);
                if(y.PE > best_single)
                {
                   best_single = y.PE;
                   best_single_index = u;
                   best_second_index = j;
                }
                ///for(int i = 0; i < num_node; i++)
                ///{
                ///    cout << y.satisfied[i];
                ///}
                ///cout << "\n";
                cout << y.PE << endl;
             }
             u++;
          }
          cout << "indices : " << best_single_index << " ";
          cout << best_second_index;
          cout << " best : " << best_single << endl;
          */





        ///new nice third
        /*
        molecule y;
        int u = 250;
         while(u < 324)
         {
           for(int k = 0; k <num_node; k++)
           {
            for(int j = 0; j < num_node; j++)
            {
                for(int i = 0; i < num_node; i++)
                {
                   if(i==u || i==j || i==k)
                   {
                      y.soln[i] = 1;
                   }
                   else
                   {
                      y.soln[i] = 0;
                   }
                }
                ///for(int i = 0; i < num_node; i++)
                ///{
                ///   cout << y.soln[i];
                ///}
                ///cout << "\n";
                satisfied_vec_func(y);
                y.PE = fitness_func(y);
                if(y.PE > best_single)
                {
                   best_single = y.PE;
                   best_single_index = u;
                   best_second_index = j;
                   best_third_index = k;
                }
                ///for(int i = 0; i < num_node; i++)
                ///{
                ///    cout << y.satisfied[i];
                ///}
                ///cout << "\n";
                cout << u << " "<< y.PE << endl;
             }
           }
            u++;
          }
          cout << "indices : " << best_single_index << " " << best_second_index << " " << best_third_index<< endl;
          cout << " best : " << best_single << endl;
          */




          /*
          int co = 1;
          molecule y;
          int u = 250;
          while(u < 324)
          {
              for(int l = 0; l < num_node; l++)
              {
                 cout << " l ";
                 for(int k = l+1; k <num_node; k++)
                 {
                     cout << " k ";
                     for(int j = k+1; j < num_node; j++)
                     {
                         cout << " j ";
                         for(int i = 0; i < num_node; i++)
                         {
                              if(i==u || i==l || i==k || i==j)
                              {
                                  y.soln[i] = 1;
                              }
                              else
                              {
                                  y.soln[i] = 0;
                              }
                         }
                         ///for(int i = 0; i < num_node; i++)
                         ///{
                         ///    cout << y.soln[i];
                         ///}
                         ///cout << "\n";
                         satisfied_vec_func(y);
                         y.PE = fitness_func(y);
                         if(y.PE > best_single)
                         {
                             best_single = y.PE;
                             best_single_index = u;
                             best_second_index = j;
                             best_third_index = k;
                             best_fourth_index = l;
                         }
                         if(best_single>=12106)
                         {
                             goto pr;
                         }
                          ///for(int i = 0; i < num_node; i++)
                         ///{
                         ///    cout << y.satisfied[i];
                         ///}
                         ///cout << "\n";
                         cout << co << "  " <<  u << " "<< y.PE << endl;
                         co++;
                      }
                  }
             }
             u++;
          }
          pr: cout << "indices : " << best_single_index << " " << best_second_index << " " << best_third_index <<" " << best_fourth_index << endl;
          cout << " best : " << best_single << endl;
          */






            /*
            ///total third
            int vc = 0;
            molecule z;
            int w = 0;
            while(w < num_node)
            {
              for(int j = 0; j < num_node; j++)
              {
                 for(int p = 0; p < num_node; p++)
                 {
                    for(int i = 0; i < num_node; i++)
                    {
                       if(i==w || i==p || i==j)
                       {
                           z.soln[i] = 1;
                       }
                       else
                       {
                          z.soln[i] = 0;
                       }
                    }
                    for(int i = 0; i < num_node; i++)
                    {
                      cout << z.soln[i];
                    }
                    cout << "\n";
                    satisfied_vec_func(z);
                    z.PE = fitness_func(z);
                    if(z.PE > best_single)
                    {
                       best_single = z.PE;
                       best_single_index = w;
                       best_second_index = p;
                       best_third_index = j;
                    }
                    cout << vc << " " << z.PE << endl;
                 }
              }
              w++;
              vc++;
            }

            cout <<"\n\n";
            cout << "index :" << best_single_index << ", " << best_second_index << ", " << best_third_index << " ";
            cout << "Result : " << best_single;
            */




            /*
            molecule h;
            for(int j = 0; j < num_node; j++)
            {
              for(int i = 0; i < num_node; i++)
              {
                  if(i==best_single_index || i==best_second_index || i==j)
                  {
                      h.soln[i] = 1;
                  }
                  else
                  {
                      h.soln[i] = 0;
                  }
              }

              satisfied_vec_func(h);
              h.PE = fitness_func(h);
              if(h.PE > value)
              {
                 value = h.PE;
                 ind = j;
              }
              cout << h.PE << endl;

            }

            cout << "index : " << ind << " ";
            cout << "Result : " << value << endl;
            */





            /*
            molecule t;
            for(int j = 0; j < num_node; j++)
            {
              for(int i = 0; i < num_node; i++)
              {
                  if(i==best_single_index || i==best_second_index || i==ind || i==j)
                  {
                      t.soln[i] = 1;
                  }
                  else
                  {
                      t.soln[i] = 0;
                  }
              }

              satisfied_vec_func(t);
              t.PE = fitness_func(t);
              if(t.PE > value1)
              {
                 value1 = t.PE;
                 ind1 = j;
              }
              cout << t.PE << endl;

            }

            cout << "index : " << ind1 << " ";
            cout << "Result : " << value1;
            */




            /*
            molecule s;
            for(int j = 0; j < num_node; j++)
            {
              for(int i = 0; i < num_node; i++)
              {
                  if(i==best_single_index || i==best_second_index || i==ind || i==ind1 || i==j)
                  {
                      s.soln[i] = 1;
                  }
                  else
                  {
                      s.soln[i] = 0;
                  }
              }

              satisfied_vec_func(s);
              s.PE = fitness_func(s);
              if(s.PE > value2)
              {
                 value2 = s.PE;
                 ind2 = j;
              }
              cout << s.PE << endl;

            }

            cout << "index : " << ind2 << " ";
            cout << "Result : " << value2;
            */




           /*
           int we = 1;
           for(int i = 0; i < 100; i++)
           {
              for(int l = 0; l < 100; l++)
              {
               for(int j = 0; j < 100; j++)
               {
                for(int k = 0; k < 100; k++)
                {
                    if(k == i || k==j|| k==l )
                    {
                        cout << "1";
                    }
                    cout << "0";
                }
                cout << "  " << we;
                we++;
                cout << "\n";
               }
              }

           }*/



            /*
           for(int i = 0; i < 500; i++)
           {
               if(i==0 || i==1 || i==2 || i==3 || i==4)
               {
                   cout << "1";
               }
               else{
                   cout <<"0";
               }
           }
           */





            int coo = 1;
            molecule x;
            istringstream iss;
            ifstream get_data; //Get distance values from problem set

		    get_data.open("data.txt");
            string line;
            while(!get_data.eof())
            {
                getline(get_data, line);
                if(get_data.good())
                {
                    //cout << line << endl;
                    iss.clear();
                    iss.str(line);
                    for(int i = 0 ; i < line.length(); i++)
                    {
                        //cout << line[i] << endl;
                        //data[i] = line[i] - '0';
                        x.soln[i] = line[i] - '0';
                    }
                    //cout << endl;
                    ///for(int i = 0; i < num_node; i++)
                    ///{
                    ///    cout << x.soln[i];
                    ///}
                    ///cout << "\n";
                    if(coo>=3465567)
                    {
                        for(int i = 0; i < num_node; i++)
                        {
                            cout << x.soln[i];
                        }
                        cout << "\n";
                    }

                    satisfied_vec_func(x);
                    x.PE = fitness_func(x);
                    cout << coo << " " << x.PE;
                    cout << endl;
                    if(x.PE > value)
                    {
                       value = x.PE;
                       for(int i = 0; i < num_node; i++)
                       {
                           data[i] = x.soln[i];
                       }
                    }
                    coo++;

                }
            }
            cout << "\n\n";
            cout << "Result = " << value << endl;

            for(int i = 0; i < num_node; i++)
            {
                if(data[i]==1)
                {
                    cout << "Index : "<< i << " ";
                }
            }








          /*
          molecule e;
          for(int j = 0; j < num_node; j++)
          {
            for(int i = 0; i < num_node; i++)
            {
                if(i==71)
                {
                    e.soln[i] = 1;
                }
                else
                {
                    e.soln[i] = 0;
                }
            }
            satisfied_vec_func(e);
            e.PE = fitness_func(e);
            if(e.PE > value)
            {
                value = e.PE;
                ind = j;
            }
          }
          cout << "\n\n";
          cout << value;
          */





       return 0;

}

