#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>


using namespace std;

//==============================================================================
//Функция рандома:
int RandomValue(int MaxValue)
{
	long int r;
	r = rand()%(MaxValue+1);
	return r;
}
//------------------------------------------------------------------------------
//Функция подсчета длинны гена:
int LengthGen(float MinValue, float MaxValue, float Accuracy)
{
	int h, Length;
	h = (MaxValue - MinValue)/(Accuracy/10.);   //Интервал разбивается на h равных частей.
	Length = log(h+1.)/log(2.)+1.;
	return Length;
}
//------------------------------------------------------------------------------
//Функция подсчета Y:
double FunctionY(int NumberOfFunction, vector <double> &XVector)
{
	double Y;
	switch (NumberOfFunction)
	{
		case 0:   		//Функция Растригина
		{
			Y = 0.1*XVector[0]*XVector[0]+0.1*XVector[1]*XVector[1]-4.0*cos(0.8*XVector[0])-4.0*cos(0.8*XVector[1])+8.0;
			break;
		}
	}
	return Y;
}
//------------------------------------------------------------------------------
//Функция подсчета Y-Yist:
bool FunctionIsSolved(int NumberOfFunction, vector <double> &XVector, float Accuracy)
{
	bool IsSolved;
	IsSolved = false;
	switch (NumberOfFunction)
	{
		case 0:   		//Функция Растригина
		{
			if ((fabs(XVector[0]-0.0) < Accuracy) && (fabs(XVector[1]-0.0) < Accuracy))
			{
				IsSolved = true;
			}
			break;
		}
	}
	return IsSolved;
}
//==============================================================================
//Класс Индивид:
class CIndivid
{
	public:
		//Переменные класса:
			int LengthValue, NumberOfVariables;
			int SType, RType, MType;				//самонастройка: тип селекции, рек-ии и мутации
			float MinValue, MaxValue, AccuracyValue;
			double YValue, FitnessValue;

		//Векторы класса:
			std::vector <bool> GenotypeVector;     //вектор генотипа, двоичное представление решения
			std::vector <double> XVector;          //вектор переменных, вещ-ое представление решения

		//Функции класса:
			CIndivid::CIndivid(float MIN, float MAX, float Accuracy, int NumOfVars, int Length);
			CIndivid::~CIndivid();
			void operator = (CIndivid A);
			void FunctionValue(int NumberOfFunction); 	//расчет иксов и игрика
			void Fitness(double YMinValue);			//расчет пригодности
			//void Operators(double SProbability, double RProbability, double MProbability);                         //самонастрока - выбор операторов
			void Mutation(double Probability);
};
//------------------------------------------------------------------------------
CIndivid::CIndivid(float MIN, float MAX, float Accuracy, int NumOfVars, int Length)
{
	LengthValue = Length;
	NumberOfVariables = NumOfVars;
	MinValue = MIN;
	MaxValue = MAX;
	AccuracyValue = Accuracy;

	for (int i = 0; i < LengthValue*NumberOfVariables; i++)
	{
		GenotypeVector.push_back(RandomValue(1));
	}
	for (int i = 0; i < NumberOfVariables; i++)
	{
		XVector.push_back(0);
	}
}
//------------------------------------------------------------------------------
CIndivid::~CIndivid()
{

}
//------------------------------------------------------------------------------
void CIndivid::operator = (CIndivid A)
{
	for (int i = 0; i < LengthValue*NumberOfVariables; i++)
	{
		GenotypeVector[i] = A.GenotypeVector[i];
	}
	for (int i = 0; i < NumberOfVariables; i++)
	{
		XVector[i] = A.XVector[i];
	}
	YValue = A.YValue;
	FitnessValue = A.FitnessValue;
	/*SType = A.SType;
	RType = A.RType;
	MType = A.MType;   */
}
//------------------------------------------------------------------------------
void CIndivid::FunctionValue(int NumberOfFunction)
{
	for (int i = 0; i < NumberOfVariables; i++)
	{
		XVector[i] = 0;
	}
	for (int i = 0; i < NumberOfVariables; i++)
	{
		for (int j = 0; j < LengthValue; j++)
		{
			XVector[i] += GenotypeVector[LengthValue*(i+1)-1-j]*pow(2.,j);
		}
		XVector[i] = MinValue + XVector[i]*(AccuracyValue/10.);
		if (XVector[i] > MaxValue)
		{
			XVector[i] = XVector[i] - fabs(MaxValue - MinValue);
		}
	}
	YValue = FunctionY(NumberOfFunction, XVector);
}
//------------------------------------------------------------------------------
void CIndivid::Fitness(double YMinValue)
{
	double YV;

	YV = YValue + YMinValue;
	if (YV == 0)
	{
		FitnessValue = 1./(YV+0.0000000001);
	}
	else FitnessValue = 1./YV;
}
//------------------------------------------------------------------------------
void CIndivid::Mutation(double Probability)
{
	double r;

	for (int i = 0; i < LengthValue*NumberOfVariables; i++)
	{
		r = RandomValue (10000);
		r = r/10000.;
		if (r < Probability)
		{
			if (GenotypeVector[i] == 1)
			{
				GenotypeVector[i] = 0;
			}
			else GenotypeVector[i] = 1;
		}
	}
}
//==============================================================================
//Класс Популяция:
class CPopulation
{
	public:
			//Переменные класса:
			int PopulationSize, NumberOfVariables, LengthValue;
			int CurrentGeneration, BestGeneration;
			int Selections, Recombinations, Mutations;
			float MinValue, MaxValue, AccuracyValue;
			double MinY;

			//Векторы класса:
			std::vector <CIndivid> PopulationVector;		//Вектор индивидов
			std::vector <CIndivid> ParentsVector;           //Вектор родителей
			std::vector <CIndivid> OffspringsPopulation;
			std::vector <double> BestFitnessVector;         //Вектор лучших значений пригодности (на каждом поколении)
			std::vector <double> AverageFitnessVector;      //Вектор средних значений пригодности (на каждом поколении)
			std::vector <double> WorstFitnessVector;        //Вектор худших значений пригодности (на кождом поколении)
			std::vector <double> SelectionProbability;		// 0 - турнирная, 1 - пропорциональная, 2 - ранговая
			std::vector <double> RecombinationProbability;  // 0 - одноточ, 1 - двухточ, 3 - равномерное
			std::vector <double> MutationProbability;       // 0 - слабая, 2 - средняя, 3 - сильная
			
			//Функции класса:
			CPopulation::CPopulation(int SizeOfPopulation, float MIN, float MAX, float Accuracy, int NumOfVars, int Length, int NumOfSs, int NumOfRs, int NumOfMs);
			CPopulation::~CPopulation();
			void Fitness(int NumberOfFunction, vector <CIndivid> &Individs);
			void StatFitness();
			void Selection();
			void Recombination();
			void Mutation();
			void ComputeProbs(int NumberOfGenerations);							//самонастройка - пересчет вер-тей
			void NewGeneration(int TypeOfForming);
};
//------------------------------------------------------------------------------
CPopulation::CPopulation(int SizeOfPopulation, float MIN, float MAX, float Accuracy, int NumOfVars, int Length, int NumOfSs, int NumOfRs, int NumOfMs)
{
	PopulationSize = SizeOfPopulation;
	MinValue = MIN;
	MaxValue = MAX;
	AccuracyValue = Accuracy;
	NumberOfVariables = NumOfVars;
	LengthValue = Length;
	Selections = NumOfSs;
	Recombinations = NumOfRs;
	Mutations = NumOfMs;

	CurrentGeneration = 0;
	BestGeneration = 1;

	for (int i = 0; i < PopulationSize+2; i++)
	{
		PopulationVector.push_back(CIndivid(MinValue, MaxValue, AccuracyValue, NumberOfVariables, LengthValue));
	}
	PopulationVector[PopulationSize].FitnessValue = 0.;			//лучший в поколении
	PopulationVector[PopulationSize+1].FitnessValue = 0.;       //лучший за все время

	double s,r,m;            									//установка равных вероятностей:
	s = 1./Selections;
	r = 1./Recombinations;
	m = 1./Mutations;
	for (int i = 0; i < Selections; i++)                        //селекции
	{
		SelectionProbability.push_back(s);
	}
	for (int i = 0; i < Recombinations; i++)                    //рекомбинации
	{
		RecombinationProbability.push_back(r);
	}
	for (int i = 0; i < Mutations; i++)                         //мутации
	{
		MutationProbability.push_back(m);
	}

	ParentsVector.clear();
	OffspringsPopulation.clear();
	BestFitnessVector.clear();
	AverageFitnessVector.clear();
	WorstFitnessVector.clear();

}
//------------------------------------------------------------------------------
CPopulation::~CPopulation()
{

}
//------------------------------------------------------------------------------
void CPopulation::Fitness(int NumberOfFunction, vector <CIndivid> &Individs)
{
	double YMinValue;
	//double TotalFitness, BestFitness, AverageFitness, WorstFitness;

	//PopulationVector[PopulationSize].FitnessValue = 0;

	for (int i = 0; i < PopulationSize; i++)
	{
		Individs[i].FunctionValue(NumberOfFunction);
	}
	YMinValue = 0;
	for (int i = 0; i < PopulationSize; i++)
	{
		if (Individs[i].YValue < YMinValue)
		{
			YMinValue = Individs[i].YValue;
		}
	}
	for (int i = 0; i < PopulationSize; i++)
	{
		Individs[i].Fitness(YMinValue);
		/*if (Individs[i].FitnessValue > Individs[PopulationSize].FitnessValue)
		{
			Individs[PopulationSize] = Individs[i];
		}*/
	}
	/*
	if (Individs[PopulationSize].FitnessValue > Individs[PopulationSize+1].FitnessValue)
	{
		Individs[PopulationSize+1] = Individs[PopulationSize];
		BestGeneration = CurrentGeneration+1;
	} */
	/*
	TotalFitness = 0.;
	for (int i = 0; i < PopulationSize; i++)
	{
		TotalFitness += PopulationVector[i].FitnessValue;		//общая пригодность
	}
	AverageFitness = TotalFitness/PopulationSize;				//средняя пригодность
	AverageFitnessVector.push_back(AverageFitness);
	BestFitness = PopulationVector[PopulationSize].FitnessValue;	//лучшая пригодность
	BestFitnessVector.push_back(BestFitness);
	WorstFitness = PopulationVector[0].FitnessValue;
	for (int i = 1; i < PopulationSize; i++)
	{
		if (PopulationVector[i].FitnessValue < WorstFitness)
		{
			WorstFitness = PopulationVector[i].FitnessValue;		//худшая пригодность
		}
	}
	WorstFitnessVector.push_back(WorstFitness);*/
}
//------------------------------------------------------------------------------
void CPopulation::StatFitness()
{
	double TotalFitness, BestFitness, AverageFitness, WorstFitness;

	PopulationVector[PopulationSize].FitnessValue = 0;
	for (int i = 0; i < PopulationSize; i++)
	{
		if (PopulationVector[i].FitnessValue > PopulationVector[PopulationSize].FitnessValue)
		{
			PopulationVector[PopulationSize] = PopulationVector[i];
		}
	}
	if (PopulationVector[PopulationSize].FitnessValue > PopulationVector[PopulationSize+1].FitnessValue)
	{
		PopulationVector[PopulationSize+1] = PopulationVector[PopulationSize];
		BestGeneration = CurrentGeneration+1;
	}
	TotalFitness = 0.;
	for (int i = 0; i < PopulationSize; i++)
	{
		TotalFitness += PopulationVector[i].FitnessValue;		//общая пригодность
	}
	AverageFitness = TotalFitness/PopulationSize;				//средняя пригодность
	AverageFitnessVector.push_back(AverageFitness);
	BestFitness = PopulationVector[PopulationSize].FitnessValue;	//лучшая пригодность
	BestFitnessVector.push_back(BestFitness);
	WorstFitness = PopulationVector[0].FitnessValue;
	for (int i = 1; i < PopulationSize; i++)
	{
		if (PopulationVector[i].FitnessValue < WorstFitness)
		{
			WorstFitness = PopulationVector[i].FitnessValue;		//худшая пригодность
		}
	}
	WorstFitnessVector.push_back(WorstFitness);

}
//------------------------------------------------------------------------------
void CPopulation::Selection()
{
	std::vector <double> SelV;
	double rns, sel;
	int TournamentSize;

	ParentsVector.clear();
	OffspringsPopulation.clear();

	sel = 0;
	for (int s = 0; s < Selections; s++)
	{
		sel += SelectionProbability[s];
		SelV.push_back(sel);
	}

	ofstream fout;
	fout.open("test.txt", ios::app);
	fout<<"Probabilities: (";
	for (int i = 0; i < SelV.size(); i++)
	{
		if (i < SelV.size()-1)
		{
			fout<<SelV[i]<<";";
		}
		else
		{
			fout<<SelV[i]<<")\n";
		}
	}
	
	for (int i = 0; i < PopulationSize; i++)
	{
		OffspringsPopulation.push_back(CIndivid(MinValue, MaxValue, AccuracyValue, NumberOfVariables, LengthValue));
		rns = RandomValue(10000);
		rns = rns/10000.;
		for (int s = 0; s < SelV.size(); s++)
		{
			if (rns <= SelV[s])
			{
				OffspringsPopulation[i].SType = s;
				s = SelV.size();
			}
		}
	}
	//fout<<"Selection for parents: ";
	for (int i = 0; i < 2*PopulationSize; i++)
	{
		ParentsVector.push_back(CIndivid(MinValue, MaxValue, AccuracyValue, NumberOfVariables, LengthValue));
	}

	std::vector <CIndivid> SortedPopulation;
	std::vector <int> IntSortedPopulation;
	CIndivid A (MinValue, MaxValue, AccuracyValue, NumberOfVariables, LengthValue);
	int IntA;
	for (int j = 0; j < PopulationSize; j++)
	{
		SortedPopulation.push_back(CIndivid(MinValue, MaxValue, AccuracyValue, NumberOfVariables, LengthValue));
		SortedPopulation[j] = PopulationVector[j];
		IntSortedPopulation.push_back(j);
	}
	for (int j = 0; j < SortedPopulation.size(); j++)
	{
		for (int k = 0; k < SortedPopulation.size()-1; k++)
		{
			if (SortedPopulation[k].FitnessValue > SortedPopulation[k+1].FitnessValue)
			{
				A = SortedPopulation[k];
				SortedPopulation[k] = SortedPopulation[k+1];
				SortedPopulation[k+1] = A;
				IntA = IntSortedPopulation[k];
				IntSortedPopulation[k] = IntSortedPopulation[k+1];
				IntSortedPopulation[k+1] = IntA;
			}
		}
	}
	/*fout<<"\n\nSorted Population:";
	for (int j = 0; j < SortedPopulation.size(); j++)
	{
		fout<<"\n"<<j<<"\t";
		for (int k = 0; k < SortedPopulation[j].GenotypeVector.size(); k++)
		{
			fout<<SortedPopulation[j].GenotypeVector[k];
		}
		fout<<"\t"<<SortedPopulation[j].FitnessValue;
	}
	*/
	for (int i = 0; i < PopulationSize; i++)
	{
		int rt = OffspringsPopulation[i].SType;
		//fout<<"\n"<<i<<"\t"<<ParentsSelection[i];
		switch (OffspringsPopulation[i].SType)
		{
			case 0:			//турнирная
			{
				//fout<<" - tournament";
				TournamentSize = 2;
				for (int s = 0; s < 2; s++)
				{
					int Number, BestNumber;
					std::vector <int> TournamentSelection;
					std::vector <int> NumbersIndivids;

					TournamentSelection.clear();

					for (int j = 0; j < PopulationSize; j++)
					{
						NumbersIndivids.push_back(j);
					}
					for (int j = 0; j < TournamentSize; j++)
					{
						Number = RandomValue(NumbersIndivids.size()-1);
						TournamentSelection.push_back(NumbersIndivids[Number]);
						NumbersIndivids.erase(NumbersIndivids.begin()+Number);
					}

					BestNumber = TournamentSelection[0];
					for (int j = 0; j < TournamentSelection.size(); j++)
					{
						if (PopulationVector[TournamentSelection[j]].FitnessValue > PopulationVector[BestNumber].FitnessValue)
						{
							BestNumber = TournamentSelection[j];
						}
					}
					/*fout<<"\n\tWinner "<<s+1<<" (";
					for (int k = 0; k < TournamentSelection.size(); k++)
					{
						if (k < TournamentSelection.size()-1)
						{
							fout<<TournamentSelection[k]<<";";
						}
						else
						{
							fout<<TournamentSelection[k]<<")";
						}
					}
					fout<<" = "<<BestNumber<<"\t";
					*/
					ParentsVector[2*i+s] = PopulationVector[BestNumber];
					/*for (int k = 0; k < ParentsVector[2*i+s].GenotypeVector.size(); k++)
					{
						fout<<ParentsVector[2*i+s].GenotypeVector[k];
					}
					*/
				}
				break;
			}
			case 1:			//турнирная
			{
				//fout<<" - tournament";
				TournamentSize = 5;
				for (int s = 0; s < 2; s++)
				{
					int Number, BestNumber;
					std::vector <int> TournamentSelection;
					std::vector <int> NumbersIndivids;

					TournamentSelection.clear();

					for (int j = 0; j < PopulationSize; j++)
					{
						NumbersIndivids.push_back(j);
					}
					for (int j = 0; j < TournamentSize; j++)
					{
						Number = RandomValue(NumbersIndivids.size()-1);
						TournamentSelection.push_back(NumbersIndivids[Number]);
						NumbersIndivids.erase(NumbersIndivids.begin()+Number);
					}

					BestNumber = TournamentSelection[0];
					for (int j = 0; j < TournamentSelection.size(); j++)
					{
						if (PopulationVector[TournamentSelection[j]].FitnessValue > PopulationVector[BestNumber].FitnessValue)
						{
							BestNumber = TournamentSelection[j];
						}
					}
					/*fout<<"\n\tWinner "<<s+1<<" (";
					for (int k = 0; k < TournamentSelection.size(); k++)
					{
						if (k < TournamentSelection.size()-1)
						{
							fout<<TournamentSelection[k]<<";";
						}
						else
						{
							fout<<TournamentSelection[k]<<")";
						}
					}
					fout<<" = "<<BestNumber<<"\t";
					*/
					ParentsVector[2*i+s] = PopulationVector[BestNumber];
					/*for (int k = 0; k < ParentsVector[2*i+s].GenotypeVector.size(); k++)
					{
						fout<<ParentsVector[2*i+s].GenotypeVector[k];
					}
					*/
				}
				break;
			}
			case 2:			//турнирная
			{
				//fout<<" - tournament";
				TournamentSize = 7;
				for (int s = 0; s < 2; s++)
				{
					int Number, BestNumber;
					std::vector <int> TournamentSelection;
					std::vector <int> NumbersIndivids;

					TournamentSelection.clear();

					for (int j = 0; j < PopulationSize; j++)
					{
						NumbersIndivids.push_back(j);
					}
					for (int j = 0; j < TournamentSize; j++)
					{
						Number = RandomValue(NumbersIndivids.size()-1);
						TournamentSelection.push_back(NumbersIndivids[Number]);
						NumbersIndivids.erase(NumbersIndivids.begin()+Number);
					}

					BestNumber = TournamentSelection[0];
					for (int j = 0; j < TournamentSelection.size(); j++)
					{
						if (PopulationVector[TournamentSelection[j]].FitnessValue > PopulationVector[BestNumber].FitnessValue)
						{
							BestNumber = TournamentSelection[j];
						}
					}
					/*fout<<"\n\tWinner "<<s+1<<" (";
					for (int k = 0; k < TournamentSelection.size(); k++)
					{
						if (k < TournamentSelection.size()-1)
						{
							fout<<TournamentSelection[k]<<";";
						}
						else
						{
							fout<<TournamentSelection[k]<<")";
						}
					}
					fout<<" = "<<BestNumber<<"\t";
					*/
					ParentsVector[2*i+s] = PopulationVector[BestNumber];
					/*for (int k = 0; k < ParentsVector[2*i+s].GenotypeVector.size(); k++)
					{
						fout<<ParentsVector[2*i+s].GenotypeVector[k];
					}
					*/
				}
				break;
			}
			case 3:			//пропорциональная
			{
				//fout<<" - proportional\t";

				std::vector <double> ProportionalSelection;
				double TotalFitness, pv, ps;

				TotalFitness = 0;
				for (int j = 0; j < PopulationSize; j++)
				{
					TotalFitness += PopulationVector[j].FitnessValue;
				}
				ps = 0;
				for (int j = 0; j < PopulationSize; j++)
				{
					pv = PopulationVector[j].FitnessValue/TotalFitness;
					ps += pv;
					ProportionalSelection.push_back(ps);
				}
				/*fout<<"(";
				for (int k = 0; k < ProportionalSelection.size(); k++)
				{
					if (k < ProportionalSelection.size()-1)
					{
						fout<<ProportionalSelection[k]<<";";
					}
					else fout<<ProportionalSelection[k]<<")";
				}
				*/
				for (int s = 0; s < 2; s++)
				{
					ps = RandomValue(10000)/10000.;
					for (int j = 0; j < ProportionalSelection.size(); j++)
					{
						if (ps < ProportionalSelection[j])
						{
							ParentsVector[2*i+s] = PopulationVector[j];
							j = ProportionalSelection.size();
						}
					}
					/*fout<<"\n\t";
					for (int k = 0; k < ParentsVector[2*i+s].GenotypeVector.size(); k++)
					{
						fout<<ParentsVector[2*i+s].GenotypeVector[k];
					}
					*/
				}
				break;
			}
			case 4:			//ранговая
			{
				//fout<<" - rank\n\t";

				std::vector <double> RankPopulation;
				std::vector <double> SortedRanks;
				std::vector <double> RankSelection;
				int num;
				double pv, ps, r, TotalRank;

				for (int j = 0; j < SortedPopulation.size(); j++)
				{
					SortedRanks.push_back(j+1);
					RankPopulation.push_back(0);
				}
				/*fout<<"pre-Ranks:\t";
				for (int j = 0; j < RankPopulation.size(); j++)
				{
					fout<<RankPopulation[j]<<" ";
				}
				fout<<"\n\tRanks:\t\t";
				*/
				for (int j = 0; j < SortedPopulation.size()-1; j++)
				{
					r = SortedRanks[j];
					num = 1;
					for (int k = j+1; k < SortedPopulation.size(); k++)
					{
						if (SortedPopulation[j].FitnessValue == SortedPopulation[k].FitnessValue)
						{
							r += SortedRanks[k];
							num++;
						}
					}
					r = r/num;
					for (int k = j; k < j+num; k++)
					{
						SortedRanks[k] = r;
					}
					j = j+num-1;
				}
				/*for (int j = 0; j < RankPopulation.size(); j++)
				{
					fout<<RankPopulation[j]<<" ";
				}
				*/
				for (int j = 0; j < SortedRanks.size(); j++)
				{
					RankPopulation[IntSortedPopulation[j]] = SortedRanks[j];
				}
				TotalRank = 0;
				for (int j = 0; j < RankPopulation.size(); j++)
				{
					TotalRank += RankPopulation[j];
				}
				pv = 0;
				for (int j = 0; j < RankPopulation.size(); j++)
				{
					pv += RankPopulation[j]/TotalRank;
					RankSelection.push_back(pv);
				}
				/*fout<<"\n\t(";
				for (int j = 0; j < RankSelection.size(); j++)
				{
					if (j < RankSelection.size()-1)
					{
						fout<<RankSelection[j]<<";";
					}
					else fout<<RankSelection[j]<<")";
				}
				*/
                for (int s = 0; s < 2; s++)
				{
					ps = RandomValue(10000)/10000.;
					for (int j = 0; j < RankSelection.size(); j++)
					{
						if (ps < RankSelection[j])
						{
							ParentsVector[2*i+s] = PopulationVector[j];
							j = RankSelection.size();
						}
					}
					/*fout<<"\n\t";
					for (int k = 0; k < ParentsVector[2*i+s].GenotypeVector.size(); k++)
					{
						fout<<ParentsVector[2*i+s].GenotypeVector[k];
					}
					*/
				}
				break;
			}
		}
	}
	/*fout<<"\n\n\tParents:\n";
	for (int i = 0; i < ParentsVector.size(); i++)
	{
		fout<<"\n"<<i+1<<"\t";
		for (int j = 0; j < ParentsVector[i].GenotypeVector.size(); j++)
		{
			fout<<ParentsVector[i].GenotypeVector[j];
		}
	}
	*/
	fout.close();
}
//------------------------------------------------------------------------------
void CPopulation::Recombination()
{
	std::vector <double> RecV;
	double rec, rnr;

	rec = 0;
	for (int i = 0; i < RecombinationProbability.size(); i++)
	{
		rec += RecombinationProbability[i];
		RecV.push_back(rec);
	}

	ofstream fout;
	fout.open("test.txt", ios::app);
	fout<<"Probabilities: (";
	for (int i = 0; i < RecV.size(); i++)
	{
		if (i < RecV.size()-1)
		{
			fout<<RecV[i]<<";";
		}
		else
		{
			fout<<RecV[i]<<")\n";
		}
	}
	for (int i = 0; i < PopulationSize; i++)
	{
		rnr = RandomValue(10000);
		rnr = rnr/10000.;
		for (int r = 0; r < RecV.size(); r++)
		{
			if (rnr < RecV[r])
			{
				OffspringsPopulation[i].RType = r;
				r = RecV.size();
			}
		}
	}
	for (int i = 0; i < ParentsVector.size(); i+=2) 
	{
		switch (OffspringsPopulation[i/2.].RType) 
		{
			case 0:			//одноточечное
			{
				int BreakPoint, rv;
				std::vector <CIndivid> Offsprings;

				for (int j = 0; j < 2; j++)
				{
					Offsprings.push_back(CIndivid(MinValue, MaxValue, AccuracyValue, NumberOfVariables, LengthValue));
				}
				BreakPoint = RandomValue(LengthValue*NumberOfVariables-2)+1;
				for (int j = 0; j < 2; j++) 
				{
					for (int k = 0; k < BreakPoint; k++) 
					{
						Offsprings[j].GenotypeVector[k] = ParentsVector[i+j].GenotypeVector[k];
					}
					for (int k = BreakPoint; k < LengthValue*NumberOfVariables; k++) 
					{
						Offsprings[j].GenotypeVector[k] = ParentsVector[i+1-j].GenotypeVector[k];
					}
				}
				rv = RandomValue(1);
				OffspringsPopulation[i/2.] = Offsprings[rv];
				/*fout<<"\nParent 1:\t";
				for (int j = 0; j < BreakPoint; j++)
				{
					fout<<ParentsVector[i].GenotypeVector[j];
				}
				fout<<"||";
				for (int j = BreakPoint; j < LengthValue*NumberOfVariables; j++)
				{
					fout<<ParentsVector[i].GenotypeVector[j];
				}
				fout<<"\nParent 2:\t";
				for (int j = 0; j < BreakPoint; j++)
				{
					fout<<ParentsVector[i+1].GenotypeVector[j];
				}
				fout<<"||";
				for (int j = BreakPoint; j < LengthValue*NumberOfVariables; j++)
				{
					fout<<ParentsVector[i+1].GenotypeVector[j];
				}
				for (int j = 0; j < Offsprings.size(); j++)
				{
					fout<<"\nChild "<<j+1<<":\t";
					for (int k = 0; k < BreakPoint; k++)
					{
						fout<<Offsprings[j].GenotypeVector[k];
					}
					fout<<"||";
					for (int k = BreakPoint; k < LengthValue*NumberOfVariables; k++) 
					{
						fout<<Offsprings[j].GenotypeVector[k];
					}
				}
				fout<<"\nOffspring:\t";
				for (int j = 0; j < OffspringsPopulation[i/2.].GenotypeVector.size(); j++) 
				{
					fout<<OffspringsPopulation[i/2.].GenotypeVector[j];
				}
				fout<<"\n";*/
				break;
			}
			case 1:         //двухточечное
			{
				int BreakPoint1, BreakPoint2, rv;
				std::vector <CIndivid> Offsprings;

				for (int j = 0; j < 2; j++)
				{
					Offsprings.push_back(CIndivid(MinValue, MaxValue, AccuracyValue, NumberOfVariables, LengthValue));
				}
				int p = LengthValue*NumberOfVariables;
				BreakPoint1 = RandomValue(LengthValue*NumberOfVariables-3)+1;
				BreakPoint2 = BreakPoint1+RandomValue(LengthValue*NumberOfVariables-BreakPoint1-2)+1;
				for (int j = 0; j < 2; j++) 
				{
					for (int k = 0; k < BreakPoint1; k++) 
					{
						Offsprings[j].GenotypeVector[k] = ParentsVector[i+j].GenotypeVector[k];
					}
					for (int k = BreakPoint1; k < BreakPoint2; k++) 
					{
						Offsprings[j].GenotypeVector[k] = ParentsVector[i+1-j].GenotypeVector[k];
					}
					for (int k = BreakPoint2; k < LengthValue*NumberOfVariables; k++) 
					{
						Offsprings[j].GenotypeVector[k] = ParentsVector[i+j].GenotypeVector[k];
					}
				}
				rv = RandomValue(1);
				OffspringsPopulation[i/2.] = Offsprings[rv]; 
				/*fout<<"\nParent 1:\t";
				for (int j = 0; j < BreakPoint1; j++) 
				{
					fout<<ParentsVector[i].GenotypeVector[j];
				}
				fout<<"||";
				for (int j = BreakPoint1; j < BreakPoint2; j++) 
				{
					fout<<ParentsVector[i].GenotypeVector[j];
				}
				fout<<"||";
				for (int j = BreakPoint2; j < LengthValue*NumberOfVariables; j++) 
				{
					fout<<ParentsVector[i].GenotypeVector[j];
				}
				fout<<"\nParent 2:\t";
				for (int j = 0; j < BreakPoint1; j++) 
				{
					fout<<ParentsVector[i+1].GenotypeVector[j];
				}
				fout<<"||";
				for (int j = BreakPoint1; j < BreakPoint2; j++) 
				{
					fout<<ParentsVector[i+1].GenotypeVector[j];
				}
				fout<<"||";
				for (int j = BreakPoint2; j < LengthValue*NumberOfVariables; j++) 
				{
					fout<<ParentsVector[i+1].GenotypeVector[j];
				}
				for (int j = 0; j < Offsprings.size(); j++) 
				{
					fout<<"\nChild "<<j+1<<":\t";
					for (int k = 0; k < BreakPoint1; k++)
					{
						fout<<Offsprings[j].GenotypeVector[k];
					}
					fout<<"||";
					for (int k = BreakPoint1; k < BreakPoint2; k++) 
					{
						fout<<Offsprings[j].GenotypeVector[k];
					}
					fout<<"||";
					for (int k = BreakPoint2; k < LengthValue*NumberOfVariables; k++) 
					{
						fout<<Offsprings[j].GenotypeVector[k];
					}
				}
				fout<<"\nOffspring:\t";
				for (int j = 0; j < OffspringsPopulation[i/2.].GenotypeVector.size(); j++) 
				{
					fout<<OffspringsPopulation[i/2.].GenotypeVector[j];
				}
				fout<<"\n";*/
				break;
			}
			case 2:         //равномерное
			{
				double rv;

				for (int j = 0; j < LengthValue*NumberOfVariables; j++)
				{
					rv = RandomValue(10000);
					rv = rv/10000.;
					if (rv < 0.5)
					{
						OffspringsPopulation[i/2.].GenotypeVector[j] = ParentsVector[i].GenotypeVector[j];
					}
					else OffspringsPopulation[i/2.].GenotypeVector[j] = ParentsVector[i+1].GenotypeVector[j];
				}
				/*fout<<"\nParent 1:\t";
				for (int j = 0; j < LengthValue*NumberOfVariables; j++)
				{
					fout<<ParentsVector[i].GenotypeVector[j];
				}
				fout<<"\nParent 2:\t";
				for (int j = 0; j < LengthValue*NumberOfVariables; j++)
				{
					fout<<ParentsVector[i+1].GenotypeVector[j];
				}
				fout<<"\nOffspring :\t";
				for (int j = 0; j < LengthValue*NumberOfVariables; j++)
				{
					fout<<OffspringsPopulation[i/2.].GenotypeVector[j];
				}
				fout<<"\n";*/
				break;
			}
		}
	}

	fout.close();
}
//------------------------------------------------------------------------------
void CPopulation::Mutation()
{
	std::vector <double> MutV;
	double mut, rnm, Probability;

	mut = 0;
	for (int i = 0; i < MutationProbability.size(); i++)
	{
		mut += MutationProbability[i];
		MutV.push_back(mut);
	}
	ofstream fout;
	fout.open("test.txt", ios::app);
	fout<<"Probabilities: (";
	for (int i = 0; i < MutV.size(); i++)
	{
		if (i < MutV.size()-1)
		{
			fout<<MutV[i]<<";";
		}
		else
		{
			fout<<MutV[i]<<")\n";
		}
	}
	for (int i = 0; i < OffspringsPopulation.size(); i++)
	{
		rnm = RandomValue(10000);
		rnm = rnm/10000.;
		for (int m = 0; m < MutV.size(); m++)
		{
			if (rnm < MutV[m])
			{
				OffspringsPopulation[i].MType = m;
				m = MutV.size();
			}
		}
	}
	for (int i = 0; i < OffspringsPopulation.size(); i++)
	{
		switch (OffspringsPopulation[i].MType)
		{
			case 0:
			{
				Probability = 1./(2.*LengthValue*NumberOfVariables);
				OffspringsPopulation[i].Mutation(Probability);
				break;
			}
			case 1:
			{
				Probability = 1./LengthValue*NumberOfVariables;
				OffspringsPopulation[i].Mutation(Probability);
				break;
			}
			case 2:
			{
				Probability = 2./LengthValue*NumberOfVariables;
				OffspringsPopulation[i].Mutation(Probability);
				break;
			}
		}
		/*fout<<"\n"<<i+1<<"\t";
		for (int j = 0; j < OffspringsPopulation[i].GenotypeVector.size(); j++)
		{
			fout<<OffspringsPopulation[i].GenotypeVector[j];
		}*/
	}
	fout.close();
}
//------------------------------------------------------------------------------
void CPopulation::ComputeProbs(int NumberOfGenerations)
{
	double PSLimit, PRLimit, PMLimit;
	double newP, AF, MaxAF, Delta;
	int Num;
	std::vector <double> OperatorsAverageFitness;

	PSLimit = 3./(10.*SelectionProbability.size());
	PRLimit = 3./(10.*RecombinationProbability.size());
	PMLimit = 3./(10.*MutationProbability.size());

	ofstream fout;
	fout.open("test.txt", ios::app);
	
	Delta = 0;
	for (int i = 0; i < SelectionProbability.size(); i++)
	{
		if (SelectionProbability[i] == PSLimit)
		{
			newP = SelectionProbability[i];
		}
		if ((SelectionProbability[i] > PSLimit) && (SelectionProbability[i] <= (PSLimit+1./(SelectionProbability.size()*NumberOfGenerations))))
		{
			newP = PSLimit;
			Delta += SelectionProbability[i]-PSLimit;
		}
		if (SelectionProbability[i] > (PSLimit+1./(SelectionProbability.size()*NumberOfGenerations)))
		{
			newP = SelectionProbability[i]-1./(SelectionProbability.size()*NumberOfGenerations);
			Delta += 1./(SelectionProbability.size()*NumberOfGenerations);
		}
		SelectionProbability[i] = newP;
	}
	OperatorsAverageFitness.clear();
	for (int i = 0; i < SelectionProbability.size(); i++) 
	{
		AF = 0;
		Num = 0;
		for (int j = 0; j < OffspringsPopulation.size(); j++)
		{
			int rt = OffspringsPopulation[j].SType;
			if (OffspringsPopulation[j].SType == i)
			{
				Num++;
				AF += OffspringsPopulation[j].FitnessValue;
			}
		}
		AF = AF/Num;
		OperatorsAverageFitness.push_back(AF);
	}
	MaxAF = OperatorsAverageFitness[0];
	Num = 0;
	for (int i = 0; i < OperatorsAverageFitness.size(); i++) 
	{
		if (OperatorsAverageFitness[i] > MaxAF) 
		{
			MaxAF = OperatorsAverageFitness[i];
			Num = i;
		}
	}
	SelectionProbability[Num] += Delta;
	fout<<"\nSelection probabilities: (";
	for (int i = 0; i < SelectionProbability.size(); i++) 
	{
		if (i < SelectionProbability.size()-1) 
		{
			fout<<SelectionProbability[i]<<";";
		}
		else fout<<SelectionProbability[i]<<")";
	}
	Delta = 0;
	for (int i = 0; i < RecombinationProbability.size(); i++)
	{
		if (RecombinationProbability[i] == PRLimit)
		{
			newP = RecombinationProbability[i];
		}
		if ((RecombinationProbability[i] > PRLimit) && (RecombinationProbability[i] < (PRLimit+1./(RecombinationProbability.size()*NumberOfGenerations))))
		{
			newP = PRLimit;
			Delta += RecombinationProbability[i]-PRLimit;
		}
		if (RecombinationProbability[i] > (PRLimit+1./(RecombinationProbability.size()*NumberOfGenerations)))
		{
			newP = RecombinationProbability[i]-1./(RecombinationProbability.size()*NumberOfGenerations);
			Delta += 1./(RecombinationProbability.size()*NumberOfGenerations);
		}
		RecombinationProbability[i] = newP;
	}
	OperatorsAverageFitness.clear();
	for (int i = 0; i < RecombinationProbability.size(); i++)
	{
		AF = 0;
		Num = 0;
		for (int j = 0; j < OffspringsPopulation.size(); j++)
		{
			int rt = OffspringsPopulation[j].RType;
			if (OffspringsPopulation[j].RType == i)
			{
				Num++;
				AF += OffspringsPopulation[j].FitnessValue;
			}
		}
		AF = AF/Num;
		OperatorsAverageFitness.push_back(AF);
	}
	MaxAF = OperatorsAverageFitness[0];
	Num = 0;
	for (int i = 0; i < OperatorsAverageFitness.size(); i++) 
	{
		if (OperatorsAverageFitness[i] > MaxAF) 
		{
			MaxAF = OperatorsAverageFitness[i];
			Num = i;
		}
	}
	RecombinationProbability[Num] += Delta;
	fout<<"\nRecombination probabilities: (";
	for (int i = 0; i < RecombinationProbability.size(); i++) 
	{
		if (i < RecombinationProbability.size()-1) 
		{
			fout<<RecombinationProbability[i]<<";";
		}
		else fout<<RecombinationProbability[i]<<")";
	}  
	Delta = 0;
	for (int i = 0; i < MutationProbability.size(); i++)
	{
		if (MutationProbability[i] == PMLimit)
		{
			newP = MutationProbability[i];
		}
		if ((MutationProbability[i] > PMLimit) && (MutationProbability[i] < (PMLimit+1./(MutationProbability.size()*NumberOfGenerations))))
		{
			newP = PMLimit;
			Delta += MutationProbability[i]-PMLimit;
		}
		if (MutationProbability[i] > (PMLimit+1./(MutationProbability.size()*NumberOfGenerations)))
		{
			newP = MutationProbability[i]-1./(MutationProbability.size()*NumberOfGenerations);
			Delta += 1./(MutationProbability.size()*NumberOfGenerations);
		}
		MutationProbability[i] = newP;
	}
	OperatorsAverageFitness.clear();
	for (int i = 0; i < MutationProbability.size(); i++) 
	{
		AF = 0;
		Num = 0;
		for (int j = 0; j < OffspringsPopulation.size(); j++) 
		{
			if (OffspringsPopulation[j].MType == i) 
			{
				Num++;
				AF += OffspringsPopulation[j].FitnessValue;
			}
		}
		AF = AF/Num;
		OperatorsAverageFitness.push_back(AF);
	}
	MaxAF = OperatorsAverageFitness[0];
	Num = 0;
	for (int i = 0; i < OperatorsAverageFitness.size(); i++) 
	{
		if (OperatorsAverageFitness[i] > MaxAF) 
		{
			MaxAF = OperatorsAverageFitness[i];
			Num = i;
		}
	}
	MutationProbability[Num] += Delta;
	fout<<"\nMutation probabilities: (";
	for (int i = 0; i < MutationProbability.size(); i++) 
	{
		if (i < MutationProbability.size()-1) 
		{
			fout<<MutationProbability[i]<<";";
		}
		else fout<<MutationProbability[i]<<")";
	}
	fout.close();
}
//------------------------------------------------------------------------------
void CPopulation::NewGeneration(int TypeOfForming)
{
	switch (TypeOfForming)
	{
		case 0:
		{
			for (int i = 0; i < PopulationSize; i++)
			{
				PopulationVector[i] = OffspringsPopulation[i];
				PopulationVector[i].SType = OffspringsPopulation[i].SType;
				PopulationVector[i].RType = OffspringsPopulation[i].RType;
				PopulationVector[i].MType = OffspringsPopulation[i].MType;
			}
			break;
		}
		case 1:
		{
			int rv;

			for (int i = 0; i < PopulationSize; i++)
			{
				PopulationVector[i] = OffspringsPopulation[i];
				PopulationVector[i].SType = OffspringsPopulation[i].SType;
				PopulationVector[i].RType = OffspringsPopulation[i].RType;
				PopulationVector[i].MType = OffspringsPopulation[i].MType;
			}
			rv = RandomValue(PopulationSize-1);
			PopulationVector[rv] = PopulationVector[PopulationSize];
			break;
		}
	}
}

int main (void)
{
	srand((unsigned)time( NULL ));
	//Переменные:
	int PopulationSize, NumberOfVariables, TournamentSize, NumberOfGenerations, NumberOfRuns;
	int TheBestGeneration;
	int Length;
	int NumS, NumR, NumM;
	int Reliability, Speed;
	float MinV, MaxV, AccuracyValue;
	bool IsSolved;
	int sv, rv, mv, ev;

	ofstream fout;
	fout.open("Statistics.txt", ios::out);
	fout.close();

	ifstream fin;
	fin.open("info.txt");
	fin>>MinV;
	fin>>MaxV;
	fin>>AccuracyValue;
	fin>>PopulationSize;
	fin>>NumberOfGenerations;
	fin>>NumberOfRuns;
	fin.close();

	cout<<"MinValue = "<<MinV;
	//cin>>MinV;
	cout<<"\tMaxValue = "<<MaxV;
	//cin>>MaxV;
	cout<<"\nAccuracy = "<<AccuracyValue;
	//cin>>AccuracyValue;
	cout<<"\nSize of population = "<<PopulationSize;
	//cin>>PopulationSize;
	cout<<"\nNumber of generations = "<<NumberOfGenerations;
	//cin>>NumberOfGenerations;
	cout<<"\tNumber of runs = "<<NumberOfRuns;
	//cin>>NumberOfRuns;
	
	NumS = 5;
	NumR = 3;
	NumM = 3;

	NumberOfVariables = 2;		//для функции растригина

	fout.open("test.txt", ios::out);
	fout.close();
	//Проверка
	Length = LengthGen (MinV, MaxV, AccuracyValue);		//вычисленние длинны вектора генотип
	CPopulation P (PopulationSize, MinV, MaxV, AccuracyValue, NumberOfVariables, Length, NumS, NumR, NumM);
	for (int g = 0; g < NumberOfGenerations; g++)
	{
		P.Fitness(0, P.PopulationVector);
		P.StatFitness();

		fout.open("test.txt", ios::app);
		fout<<"\n\n\t\tGeneration "<<g+1<<"!\n";
		cout<<"\nGeneration "<<g+1<<"!";
		/*fout<<"\n\tPopulation:\n\n";
		for (int i = 0; i < PopulationSize; i++)
		{
			fout<<i+1<<"\t";
			for (int j = 0; j < P.PopulationVector[i].GenotypeVector.size(); j++)
			{
				fout<<P.PopulationVector[i].GenotypeVector[j];
			}
			fout<<"\tFunction(";
			for (int j = 0; j < P.PopulationVector[i].XVector.size(); j++)
			{
				if (j < P.PopulationVector[i].XVector.size()-1)
				{
					fout<<P.PopulationVector[i].XVector[j]<<"; ";
				}
				else
				{
					fout<<P.PopulationVector[i].XVector.size()<<") = ";
				}
			}
			fout<<P.PopulationVector[i].YValue<<"\tFitness = "<<P.PopulationVector[i].FitnessValue<<"\n";
		}*/
		fout<<"\nBest in "<<g+1<<" generation: "<<P.PopulationVector[PopulationSize].YValue<<"\tFitness = "<<P.PopulationVector[PopulationSize].FitnessValue;
		fout<<"\nThe Best: "<<P.PopulationVector[PopulationSize+1].YValue<<"\tFitness = "<<P.PopulationVector[PopulationSize+1].FitnessValue;

		fout<<"\nSelection. ";
		fout.close();
		P.Selection();
		fout.open("test.txt", ios::app);
		fout<<"Recombination. ";
		fout.close();
		P.Recombination();
		fout.open("test.txt", ios::app);
		fout<<"Mutation. ";
		fout.close();
		P.Mutation();
		P.Fitness(0, P.OffspringsPopulation);
		/*fout.open("test.txt", ios::app);
		fout<<"\n\n\tOffsprings:\n\n";
		for (int i = 0; i < PopulationSize; i++)
		{
			fout<<i+1<<"\t";
			for (int j = 0; j < P.OffspringsPopulation[i].GenotypeVector.size(); j++)
			{
				fout<<P.OffspringsPopulation[i].GenotypeVector[j];
			}
			fout<<"\t(sel,rec,mut)=("<<P.OffspringsPopulation[i].SType<<","<<P.OffspringsPopulation[i].RType<<","<<P.OffspringsPopulation[i].MType<<")";
			fout<<"\tFunction(";
			for (int j = 0; j < P.OffspringsPopulation[i].XVector.size(); j++)
			{
				if (j < P.OffspringsPopulation[i].XVector.size()-1)
				{
					fout<<P.OffspringsPopulation[i].XVector[j]<<"; ";
				}
				else
				{
					fout<<P.OffspringsPopulation[i].XVector.size()<<") = ";
				}
			}
			fout<<P.OffspringsPopulation[i].YValue<<"\tFitness = "<<P.OffspringsPopulation[i].FitnessValue<<"\n";
		}
		fout.close();*/
		P.ComputeProbs(NumberOfGenerations);
		P.NewGeneration(0);
		/*fout.open("test.txt", ios::app);
		fout<<"\n\n\t New population:\n\n";
		for (int i = 0; i < PopulationSize; i++)
		{
			fout<<i+1<<"\t";
			for (int j = 0; j < P.PopulationVector[i].GenotypeVector.size(); j++)
			{
				fout<<P.PopulationVector[i].GenotypeVector[j];
			}
			fout<<"\n";
		}
		fout.close();*/
	}



	system("PAUSE");
	return 0;
}