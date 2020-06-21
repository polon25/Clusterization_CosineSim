//============================================================================
// Name        : main.cpp
// Author      : Jacek Pi³ka
// Description : DAMI Project -> clusterization using cosine measure
// Arguments: filePath, number of header rows, delimiter, k, filePathOut
//============================================================================

#include "data.h"
#include "cosine.h"
#include <algorithm> //For sort()
#include <list>
#include <chrono> //For timing
#include <ctime>

using namespace std;

int main (int argc, char *arg[]){
	//Start counting
	auto start = std::chrono::system_clock::now();

	//Arguments
	string filePath=arg[1];
	int headerRows=atoi(arg[2]); //Number of header lines
	char delimiter=*arg[3]; //Delimiter char
	int k=atoi(arg[4]); //k neighbors
	string filePathOut=arg[5];

	/**
	 * Importing data from file
	 */

	cout<<"Parameters:"<<endl;
	cout<<"Header's rows: "<<headerRows<<endl;
	cout<<"k:"<<k<<endl<<endl;
	cout<<"Reading file: "<<filePath<<endl;
	vector<vector<string>> dataRaw=readData(filePath,headerRows,delimiter);
	vector<string> row=dataRaw[0];
	int attributeSize=row.size();
	int dataNum=dataRaw.size();
	cout<<"Data read properly"<<endl;
	cout<<"Objects in data set: "<<dataNum<<endl<<endl;

	int n=0; //For the sake of sorting

	/**
	 * Converting vector of data to the float array
	 * with additional columns for length and cluster ID
	 */

	vector<float*> data;
	vector<float*> dataOld;

	for(int i=0; i<dataNum; i++){
		data.push_back(new float[attributeSize+3]);
		dataOld.push_back(new float[attributeSize+3]);
		float dataVector[attributeSize];
		for(int ii=0; ii<attributeSize; ii++){
			dataVector[ii]=stof(dataRaw[i][ii]);
			data[i][ii]=dataVector[ii];
			dataOld[i][ii]=dataVector[ii];
		}
		data[i][attributeSize]=vectorLength(dataVector,attributeSize);
		data[i][attributeSize+1]=-1; //Initialize claster ID value with -1
		data[i][attributeSize+2]=i; //Last cell -> original id
		dataOld[i][attributeSize]=vectorLength(dataVector,attributeSize);
		dataOld[i][attributeSize+1]=-1; //Initialize claster ID value with -1
		dataOld[i][attributeSize+2]=i; //Last cell -> original id
	}

	cout<<"Preparing data, normalization and sorting"<<endl;
	//Normalize vector and calculate angle from vector [1,0,0,...,0]
	for(int i=0; i<dataNum; i++){
		for(int ii=0; ii<attributeSize+1; ii++){
			data[i][ii]=data[i][ii]/data[i][attributeSize];
		}

		//Calculate distance from arbitrary choosen point (1,0,0,...)
		float sum=0;
		sum=(data[i][0]-1)*(data[i][0]-1);
		for(int ii=1; ii<attributeSize; ii++){
			sum+=data[i][ii]*data[i][ii];
		}
		data[i][attributeSize]=sqrt(sum);
	}
	//At the end -> 2D array of vectors + angle from [1,0,0,...,0] + claster ID + original data ID

	/**
	 * Sorting data by angle from [1,0,0,...,0] vector
	 */

	//Sorting
	vector<dataVec> dataArray;
	for(int i=0; i<dataNum;i++){
		dataArray.push_back({data[i],data[i][attributeSize]});
	}
	n = sizeof(dataArray)/sizeof(dataArray[0]);
	sort(dataArray.begin(), dataArray.end(), compareDistance);

	//Writing sorted data to tmp array
	float dataTmpArray[dataNum][attributeSize+3];
	for(int i=0; i<dataNum; i++){
		for(int ii=0; ii<attributeSize+3; ii++){
			dataTmpArray[i][ii]=dataArray[i].vector[ii];
		}
	}
	//Copying data to the "main" array
	for(int i=0; i<dataNum; i++){
		for(int ii=0; ii<attributeSize+3; ii++){
			data[i][ii]=dataTmpArray[i][ii];
		}
	}
	//At the end -> sorted 2D data float array

	auto end = std::chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end-start;
	cout<<"Elapsed time: " << 1000*elapsed_seconds.count() << "ms\n"<<endl;

	cout<<"Create neighbor matrix"<<endl;

	//Matrix of neighbors (0 - not neighbor, 1 - neighbor)
	//First id - neighbors of i-th object, Second id - objects which have i-th element as their neighbor
	typedef vector<int> row_vector;
	typedef vector<row_vector> matrix_vector;
	matrix_vector neighborMatrix(dataNum, row_vector(dataNum));

	for(int i=0; i<dataNum; i++){
		for(int ii=0; ii<dataNum; ii++){
			neighborMatrix[i][ii]=0;
		}
	}

	vector<int> objectClassTable;
	for (int i=0; i<dataNum; i++){
		objectClassTable.push_back(NOISE);
	}

	/**
	 * k+ NBC clusterization using cosine measure
	 */

	//Main algorithm's loop
	//First main loop -> searching for neighbors and making neighborMatrix
	for (int i=0; i<dataNum; i++){
		//Calculate boundary values for potential close vectors

		/**
		 * Analyzing cosine values of potential close vectors
		 * Two for loops to iterate from ith vector in both directions
		 */

		float maxDistance=0;
		int neighborCount=0;

		vector<objectInfo> closeVectors; //IDs of close vectors
		for (int ii=1; ii<max(i,dataNum-i); ii++){ //Iteration through points from the i-th object
			bool operationDone=false;//Was at least one point added to the neighborhood?
			if(i+ii<dataNum){ //Iteration 'above' i-th
				int id=i+ii;
				//Calculate euclidean distance
				float distance=0;
				for (int iii=0;iii<attributeSize;iii++){
					float a=data[i][iii]-data[id][iii];
					distance+=a*a;
				}
				distance=sqrt(distance);
				if(neighborCount<k){ //If not enough neighbors, add to the neighborhood
					objectInfo object; object.id=id; object.euclideanDistance=distance;
					closeVectors.push_back(object);
					if (distance>maxDistance){//Set the max distance
						maxDistance=distance;
					}
					neighborCount++;
					operationDone=true;
				}
				else if(distance<(maxDistance+0.0001)){ //If enough neighbors, add only if distance is smaller than max
					objectInfo object; object.id=id; object.euclideanDistance=distance;
					closeVectors.push_back(object);
					neighborCount++;
					operationDone=true;
				}
			}
			if (i-ii>=0){//Iteration 'below' i-th
				int id=i-ii;
				float distance=0;
				for (int iii=0;iii<attributeSize;iii++){
					float a=data[i][iii]-data[id][iii];
					distance+=a*a;
				}
				distance=sqrt(distance);
				if(neighborCount<k){ //If not enough neighbors, add to the neighborhood
					objectInfo object; object.id=id; object.euclideanDistance=distance;
					closeVectors.push_back(object);
					if (distance>maxDistance){//Set the max distance
						maxDistance=distance;
					}
					neighborCount++;
					operationDone=true;
				}
				else if(distance<(maxDistance+0.0001)){ //If enough neighbors, add only if distance is smaller than max
					objectInfo object; object.id=id; object.euclideanDistance=distance;
					closeVectors.push_back(object);
					neighborCount++;
					operationDone=true;
				}
			}
			if (!operationDone)//If no object was added to the neighborhood
				break;
		}

		//Sort neighboors by their cosine value
		objectInfo idAndCosine[closeVectors.size()]; //Object + cluster ID (as length, cause there's no need to create new class and functions)
		for (int ii=0;ii<(int)closeVectors.size();ii++){
			idAndCosine[ii].id=closeVectors[ii].id;
			idAndCosine[ii].euclideanDistance=closeVectors[ii].euclideanDistance;
		}
		int ni = sizeof(idAndCosine)/sizeof(idAndCosine[0]);
		sort(idAndCosine, idAndCosine+ni, compareObjects);

		//If less neighbors then k
		int kTmp=(k<ni)? k : ni;

		//Add k+ neighbors to neighborMatrix
		for(int ii=0; ii<ni; ii++){
			if (idAndCosine[ii].euclideanDistance<=(idAndCosine[kTmp-1].euclideanDistance+0.000001)){ //Check if in k+ neigborhood
				neighborMatrix[i][idAndCosine[ii].id]=1; //Add to neighborMatrix
			}
			else
				break;
		}
		//Comentary: Why I add a very small number to 'border' euclidean distance?
		//Because during calculations, two points in the exactly same spot can result with very, very small but non zero euclidean distance!
		//That's why we need a extra space for such problems
	}

	cout<<"Neighbor matrix completed"<<endl;
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	cout<<"Elapsed time: " << 1000*elapsed_seconds.count() << "ms\n";
	cout<<endl;

	//Calculate ndf and check the object class
	for (int i=0;i<dataNum;i++){
		double kn=0; double rkn=0;//k-neighoors and reversed k-neighoors

		for (int ii=0; ii<dataNum; ii++){
			kn+=neighborMatrix[i][ii];
			rkn+=neighborMatrix[ii][i];
		}

		//Check the class of object
		//Only cores can create clusters (but borders can be added to them)

		double ndf=(rkn>0)? rkn/kn : 0;

		if ((kn>=k) & (ndf>=1))
			objectClassTable[i]=CORE;
		else if ((kn>=k) & (ndf<1))
			objectClassTable[i]=BORDER;
		else
			objectClassTable[i]=NOISE;
	}

	int clasterCount=0; //"When it was clusterized" counter
	int currentClusterID=0;
	cout<<"Starting clusterization"<<endl;
	//Second main loop -> clustering objects
	//Only core points can initialize cluster, border points can be only added to already existing one
	for (int i=1; i<dataNum; i++){ //iterate through idLoop2
		//If not core or already in cluster, move to next
		if ((objectClassTable[i]!=CORE)|(data[i][attributeSize+1]>=0))
			continue;

		/**
		 * Clusterization
		 */

		list<int> corePoints;

		data[i][attributeSize+1]=currentClusterID;
		data[i][attributeSize]=clasterCount;
		clasterCount++;
		for (int ii=0; ii<dataNum; ii++){
			if (neighborMatrix[i][ii]>0){
				data[ii][attributeSize+1]=currentClusterID;
				data[ii][attributeSize]=clasterCount;
				clasterCount++;
				if (objectClassTable[ii]==CORE)//if core point
					corePoints.push_back(ii);
			}
		}

		while(corePoints.front()){//iterate through all core points in cluster
			//Check it's neighbors
			for (int iii=0;iii<dataNum;iii++){
				//If neighbor and outside cluster
				if ((neighborMatrix[corePoints.front()][iii]>0)&(data[iii][attributeSize+1]<0)){
					data[iii][attributeSize+1]=currentClusterID;
					data[iii][attributeSize]=clasterCount;
					clasterCount++;
					if (objectClassTable[iii]==CORE)//If core, add to the loop iterations
						corePoints.push_back(iii);
				}
			}
			corePoints.pop_front();
		}
		currentClusterID++;
	}

	cout<<"Clasterization Complete"<<endl;
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	cout<<"Elapsed time: " << 1000*elapsed_seconds.count() << "ms\n";
	cout<<endl;

	/*************************
	Preparing data for saving:
	Sorting by cluster
	Changing cluster ids to from 0
	*************************/

	//Sort data by cluster
	dataVec idAndCluster[dataNum]; //Object + cluster ID (as length, cause there's no need to create new class and functions)
	for (int i=0;i<dataNum;i++){
		idAndCluster[i].vector=data[i];
		idAndCluster[i].length=data[i][attributeSize+1];
	}
	n = sizeof(idAndCluster)/sizeof(idAndCluster[0]);
	sort(idAndCluster,idAndCluster+n,compareDistance);

	//Writing sorted data to tmp array
	vector<float*> dataTmp2Array; //Cause it's to big array to do it the "normal way"
	for(int i=0; i<dataNum; i++){
		dataTmp2Array.push_back(new float[dataNum]);
		for(int ii=0; ii<attributeSize+3; ii++){
			dataTmp2Array[i][ii]=idAndCluster[i].vector[ii];
		}
	}
	int currentOldID=-1;
	int currentNewID=-1;
	//Copying data to the "main" array
	//Plus changing the cluster ids (0,1,2,...)
	for(int i=0; i<dataNum; i++){
		if (currentOldID!=dataTmp2Array[i][attributeSize+1]){
			currentOldID=dataTmp2Array[i][attributeSize+1];
			currentNewID++;
		}
		for(int ii=0; ii<attributeSize+3; ii++){
			data[i][ii]=dataTmp2Array[i][ii];
		}
		data[i][attributeSize+1]=currentNewID;
	}
	cout<<"Founded clusters: "<<data[dataNum-1][attributeSize+1]+1<<endl<<endl;

	/*************************
	Saving data to the new csv file
	*************************/

	float **dataTmpSave=new float *[dataNum];
	for(int i = 0; i<dataNum; i++){
		dataTmpSave[i] = new float[attributeSize+2];
	    for (int ii=0; ii<attributeSize+1; ii++){
	    	dataTmpSave[i][ii]=dataOld[(int) data[i][attributeSize+2]][ii];
	    }
	    dataTmpSave[i][attributeSize]=data[i][attributeSize];
	    dataTmpSave[i][attributeSize+1]=data[i][attributeSize+1];
	}
	writeData(dataTmpSave,dataNum,attributeSize+2,filePathOut);
	cout<<"Data Saved"<<endl;

	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;

	cout<<"Elapsed time [all]: " << 1000*elapsed_seconds.count() << "ms\n";

	return 0;
}
