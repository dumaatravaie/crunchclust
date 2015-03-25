/*
 *  Crunch_Cluster Version 43.cpp
 *  Version: 43
 *   
 *  Created on: Dec 2, 2009
 *  Updated on: Feb 28, 2011
 *  Updated on: Mar 21, 2011
 *  Last Updated on: November 24, 2011
 *  Last Updated on: February 22, 2012
 *      Author: Dipankar
 *
 *  Please mention the name of CrunchClust in the Article if you use CrunchClust for the experiment.
 *
 *  Permission is granted for anyone to copy, use, or modify these program and document for purposes of research or education, provided this copyright notice is retained, and note is made of any changes that have been made.
 *
 *
 *  This program and the document are distributed without any warranty. All use of these programs is entirely at the user's own risk.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

using namespace std;

bool Distanceone(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

bool Distance(vector<char>& s1, vector<char>& s2, int s1_length, int s2_length,
		int seed_id_size, int test_id_size, int maxOffset, int nbdf);

bool Distancemain(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

bool Distanceonepm(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

bool Distancepm(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int seed_id_size, int test_id_size, int maxOffset,
		int nbdf);

bool Distancemainpm(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

bool Distanceonehl(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

bool Distancehl(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int seed_id_size, int test_id_size, int maxOffset,
		int nbdf);

bool Distancemainhl(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

bool Distanceonehh(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

bool Distancehh(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int seed_id_size, int test_id_size, int maxOffset,
		int nbdf);

bool Distancemainhh(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb);

int readfasta_file(ifstream &in1);

int read_file_in_memory(ifstream &in1, int &SEQ_no, int *SEQ_len,
		char *SEQ_tag[], char *SEQ_final[]);

void clean(char *seq);

void cluster();

void clusterloop();

int globalflag;

//character frequency
map<int, vector<short> > globala;

map<int, vector<short> > globalb;

//Sequence files
map<int, vector<char> > global;

map<int, vector<char> > globall;

struct seqdet {
	int length;
	int dy;
	int seq_id_size;

};

struct seqsub {
	int dy;
	int seq_id_size;
	//int clusternumber;
	int distance;

};

map<int, seqdet> dip;

map<int, seqdet> diptest;

map<int, map<int, seqsub> > dipsub;

int nbdf;

char output[200];

int keep_n = 0; // whether to retain the sequences containing N or Not ?

int n_flag = 0; // whether there is N or Not in the sequence ?
//char input[200];

int *frequency_list; //Number of sequences dereplicated for each sorted sequence.
//Or the abandance of each sequence after strict derelication

int strict = 0; //**** Whether to do strict dereplication **********************
//int strict_flag = 0;

//=============== Clustering at difference threshold ===========================
int kmin = 0;

int kmax = 0;

int ksteps = 0;

int kflag = 0;
//==============================================================================
char input_referencefile[200];

char input_sequencefile[200];

char option[20];

int test_id;

int seed_id;

int rangeed = 0;

int argflag;
//===============================================================================
void help();

char endgapflagg[10];

int endgapflag;
//===============================================================================
int loop_upper_limit;

char looplimit[10];

//int arg_size;

int loop_flag = 0;

//===============================================================================
int seqnumber;
//===============================================================================

int main(int argc, char* argv[]) {
	map<int, map<int, seqsub> >::iterator it;

	FILE* pFile;
	FILE* dereplication_name_file;
	seqdet tmp;
	seqsub tmpp;
	char ch;
	long ik = 0;
	vector<char> cd;
	//****************************  Variables For Dereplication Part ***************************************
	//=====================================================================================
	//***************** For Minimum And Maximum Sequence Length ***********************
	int minimum = 5;
	int maximum = 400;
	//*********************************************************************************

	//***************** Which means cluster the sequences which are exactly equal ( In length + 0 difference ) in length
	// ***** and then sort them according to their abundance **********************

	//=====================================================================================
	//argflag = argc;

	if (argc < 9) {
		help();
	}

	//=========================================================================================

	for (int i = 1; i < argc; i++) {

		//cout<<argv[i]<<endl;

		if (strcmp(argv[i], "--diff") == 0)
			nbdf = atoi(argv[++i]);
		else if (strcmp(argv[i], "--ref") == 0)
			strcpy(input_referencefile, argv[++i]);
		else if (strcmp(argv[i], "--in") == 0)
			strcpy(input_sequencefile, argv[++i]);
		else if (strcmp(argv[i], "--out") == 0)
			strcpy(output, argv[++i]);
		else if (strcmp(argv[i], "--d_hl") == 0)
			strcpy(option, argv[i]);
		else if (strcmp(argv[i], "--d_all") == 0)
			strcpy(option, argv[i]);
		else if (strcmp(argv[i], "--d_hl_loop") == 0)
			strcpy(option, argv[i]);
		else if (strcmp(argv[i], "--strict") == 0)
			strict = 1;
		else if (strcmp(argv[i], "--kmin") == 0) {
			kmin = atoi(argv[++i]);
			kflag = 1;
		} else if (strcmp(argv[i], "--kmax") == 0) {
			kmax = atoi(argv[++i]);
			kflag = 1;
		} else if (strcmp(argv[i], "--ksteps") == 0) {
			ksteps = atoi(argv[++i]);
			kflag = 1;
		} else if (strcmp(argv[i], "--min") == 0)
			minimum = atoi(argv[++i]);
		else if (strcmp(argv[i], "--max") == 0)
			maximum = atoi(argv[++i]);
		else if (strcmp(argv[i], "--keep_n") == 0)
			keep_n = 1;
		else if (strcmp(argv[i], "--loop") == 0) {
			loop_upper_limit = atoi(argv[++i]);
			loop_flag = 1;
		} else if (strcmp(argv[i], "--noendgaps") == 0)
			endgapflag = 0;
		else if (strcmp(argv[i], "--endgaps") == 0)
			endgapflag = 1;
		else
			help();
	}

	//cout<<argc<<endl;

	//cout<<"value of endgapflag "<<endgapflag<<endl;

	if ((argc == 9) || (strict == 1) || (kflag == 1))
		argflag = 6;
	else if (argc >= 11)
		argflag = 7;
	else
		help();

	//cout<<argc<<endl;

	//===========================================================================================

	// else {

	if (argflag == 7)

	//==================================================================================
			{
		printf(" ---- ------------ ---------- ---- ");
		printf("\n ---- CrunchClust Version 43 ---- \n");
		printf(" ---- ------------ ---------- ---- \n\n");
		printf("Elapsed time: %ld secs.\n", clock() / CLOCKS_PER_SEC);

		//====================================================================================
		//nbdf = atoi(argv[1]);
		//====================================================================================

		//strcpy(input, argv[3]);

		//strcpy(output, argv[4]);
		//strcpy(option, argv[5]);
		//========================================================================================
		/*
		 tsize = strlen(option);

		 for (int i = 12; i <= tsize; i++) {
		 looplimit[i - 12] = option[i];
		 }
		 loop_upper_limit = atoi(looplimit);
		 */

		//cout<<"hello world "<<loop_upper_limit<<"  "<<tsize<<endl;
		//char a[20] = a[];
		//found = option.find("#");
		//if (found==string::npos)
		//  help();
		// position of "live" in str
		//strcpy(looplimit, str.substr(found));
		//loop_upper_limit = atoi(looplimit);*/
		//======================================================================
		/*	strcpy(endgapflagg, argv[6]);
		 if (strcmp(endgapflagg, "--noendgaps") == 0)
		 endgapflag = 0;
		 else if (strcmp(endgapflagg, "--endgaps") == 0)
		 endgapflag = 1;
		 else {
		 printf("Please enter the correct option for endgaps \n");
		 printf(" --noendgaps    :For not counting the end gaps  \n");
		 printf(" --endgaps      :For counting the end gaps   \n");
		 exit(56);

		 }*/

		//============================================================
		int aA, cC, gG, tT;

		map<short, short> a1;

		pFile = fopen(input_referencefile, "r");

		if (pFile == NULL) {
			cout << " Input_Reference_file not found .........   " << endl;
			help();
			//exit(24);
		}

		int count = 0;
		while (!feof(pFile)) {

			//int count = 0;
			do {
				fscanf(pFile, "%c", &ch);
				if ((ch != '>') && (ik == 0) && (count == 0)) {
					cout
							<< " Your input sequence fasta file is either empty or not formatted properly . . "
							<< endl;
					cout
							<< " Please check your fasta file whether it is empty or has some blank line in it . . "
							<< endl;
					cout << " Exiting the program " << endl;
					exit(1);
				}
				if ((ch != 10) && (!isspace(ch))) {
					count++;
					cd.push_back(ch);

				}
			} while (ch != 10);

			if (!cd.empty()) {
				ik++;

				tmp.seq_id_size = count;
				//tmpp.seq_id_size = count;

			}

			count = 0;

			aA = 0;
			cC = 0;
			gG = 0;
			tT = 0;

			do {
				fscanf(pFile, "%c", &ch);
				ch = toupper(ch);

				if (ch == 'A') {
					count++;
					cd.push_back(ch);

					aA++;

				} else if (ch == 'C') {
					count++;
					cd.push_back(ch);

					cC++;

				} else if (ch == 'G') {
					count++;
					cd.push_back(ch);

					gG++;

				} else if (ch == 'T') {
					count++;
					cd.push_back(ch);

					tT++;

				}

			} while ((!feof(pFile)) && (ch != '>'));

			if (!cd.empty()) {
				//ik++;
				globala[ik].push_back(aA);
				globala[ik].push_back(cC);
				globala[ik].push_back(gG);
				globala[ik].push_back(tT);

				global[ik] = cd;

				tmp.dy = 1;
				tmp.length = count;

				//tmpp.dy = 1;

				//dipsub[count][ik] = tmpp;

				dip[ik] = tmp;
			}
			cd.clear();
			if (ch == '>') {
				cd.push_back(ch);
				count = 1;
				//break;
			}

		}

		fclose(pFile); // ferme le flux et

		//================================================================================

		//int aA, cC, gG, tT;

		//map<short, short> a1;

		ik = 0;

		pFile = fopen(input_sequencefile, "r");

		if (pFile == NULL) {
			cout << " Input Sequence file not found .........   " << endl;
			help();
			//exit(25);
		}

		count = 0;
		while (!feof(pFile)) {

			//int count = 0;
			do {
				fscanf(pFile, "%c", &ch);
				if ((ch != '>') && (ik == 0) && (count == 0)) {
					cout
							<< " Your input sequence fasta file is either empty or not formatted properly . . "
							<< endl;
					cout
							<< " Please check your fasta file whether it is empty or has some blank line in it . . "
							<< endl;
					cout << " Exiting the program " << endl;
					exit(1);
				}
				if ((ch != 10) && (!isspace(ch))) {
					count++;
					cd.push_back(ch);

				}
			} while (ch != 10);

			if (!cd.empty()) {
				ik++;

				//tmp.seq_id_size = count;
				tmpp.seq_id_size = count;

			}

			count = 0;

			aA = 0;
			cC = 0;
			gG = 0;
			tT = 0;

			do {
				fscanf(pFile, "%c", &ch);
				ch = toupper(ch);

				if (ch == 'A') {
					count++;
					cd.push_back(ch);

					aA++;

				} else if (ch == 'C') {
					count++;
					cd.push_back(ch);

					cC++;

				} else if (ch == 'G') {
					count++;
					cd.push_back(ch);

					gG++;

				} else if (ch == 'T') {
					count++;
					cd.push_back(ch);

					tT++;

				}

			} while ((!feof(pFile)) && (ch != '>'));

			if (!cd.empty()) {
				//ik++;
				globalb[ik].push_back(aA);
				globalb[ik].push_back(cC);
				globalb[ik].push_back(gG);
				globalb[ik].push_back(tT);

				globall[ik] = cd;

				//tmp.dy = 1;
				//tmp.length = count;

				tmpp.dy = 1;

				dipsub[count][ik] = tmpp;

				//dip[ik] = tmp;
			}
			cd.clear();
			if (ch == '>') {
				cd.push_back(ch);
				count = 1;
				//break;
			}

		}
		seqnumber = ik;

		fclose(pFile); // ferme le flux et

		//================================================================================
		//if (strcmp(option, "--d_hl_loop") == 0)

		if (loop_flag == 1)
			clusterloop();
		else
			cluster();

		printf("Elapsed time: %ld secs.\n", clock() / CLOCKS_PER_SEC);

	}
	//======================================================================
	else if (argflag == 6) {
		printf(" \n ---- ------------ ---------- ---- ");
		printf("\n ---- CrunchClust Version 43 ---- \n");
		printf(" ---- ------------ ---------- ---- \n\n");
		printf("Elapsed time: %ld secs.\n\n", clock() / CLOCKS_PER_SEC);

		if (strict == 1) {

			//================================================================================
			int *SEQ_len;
			int *SEQ_flag;

			char *(*SEQ_tag);
			char *(*SEQ_final);
			int Total_sequence;

			//=================================================================================
			ifstream input(input_sequencefile);
			if (!input) {
				printf("Input fasta file not found \n");
				help();
			}

			//cout << " Bingo " << Total_sequence << endl;
			Total_sequence = readfasta_file(input);
			if (Total_sequence == 0) {
				cout
						<< " Your input sequence fasta file is either empty or not formatted properly . . "
						<< endl;
				cout
						<< " Please check your fasta file whether it is empty or has some blank line in it . . "
						<< endl;
				cout << " Exiting the program .. ." << endl << endl;
				exit(1);
			}
			input.close();
			printf(
					"\nTotal number of sequences before Strict Dereplication %d \n",
					Total_sequence);

			//====================================================================
			char buffer[50];
			char output_tmp[200];
			//strcpy(output_tmp,input_sequencefile);
			strncpy(output_tmp, input_sequencefile,
					(strlen(input_sequencefile) - 6));
			output_tmp[strlen(input_sequencefile) - 6] = '\0';
			//cout<<" ----------------- "<<output_tmp<<endl;
			//sprintf(buffer, "_%d", nbdf);
			strcat(output_tmp, "_Dereplication.names");
			//cout << " output_tmp " << output_tmp << endl;
			dereplication_name_file = fopen(output_tmp, "w");
			if (dereplication_name_file == NULL) {
				cout << "Error in opening " << output_tmp << "file " << endl;
				help();
				exit(1);
			}

			//====================================================================

			if ((SEQ_len = new int[Total_sequence]) == NULL) {
				printf("Memory Problem \n");
				help();
			}
			if ((SEQ_flag = new int[Total_sequence]) == NULL) {
				printf("Memory Problem \n");
				help();
			}

			/*if ((count = new int[Total_sequence]) == NULL) {
			 printf("Memory Problem \n");
			 help();
			 } */
			//=====================================================================
			if ((SEQ_tag = new char *[Total_sequence]) == NULL) {
				printf("Memory Problem \n");
				help();
			}
			if ((SEQ_final = new char *[Total_sequence]) == NULL) {
				printf("Memory Problem \n");
				help();
			}

			//=====================================================================

			ifstream in1(input_sequencefile);
			read_file_in_memory(in1, Total_sequence, SEQ_len, SEQ_tag,
					SEQ_final);
			in1.close();

			//=============================== Strict Dereplication Begins ===========================================
			for (int i = 0; i < Total_sequence; i++) {
				SEQ_flag[i] = 0;
			}

			multimap<int, int> sort_list;
			multimap<int, int>::reverse_iterator rmi;

			for (int i = 0; i < Total_sequence; i++) {
				if (SEQ_len[i] < minimum)
					continue;

				vector<int> temp_vector;

				int counter = 1;
				int ref_flag_n = 1;

				if (SEQ_flag[i] != 1) {

					int len_ref = SEQ_len[i];

					if (len_ref > maximum) {
						len_ref = maximum;
						SEQ_len[i] = maximum;
					}

					for (int k = i + 1; k < Total_sequence; k++) {

						if (SEQ_flag[k] != 1) {
							int len_test = SEQ_len[k];
							if (len_test > maximum) {
								len_test = maximum;
								SEQ_len[k] = maximum;
							}
							if (len_test < minimum) {
								SEQ_flag[k] = 1;
								continue;
							}

							if (len_ref == len_test) {
								int test_flag = 1;

								for (int p = 0; p < len_ref; p++) {

									if (SEQ_final[i][p] != SEQ_final[k][p]) {
										test_flag = 0;
										break;
									}

								}
								if (test_flag == 1) {
									counter++;
									SEQ_flag[k] = 1;
									temp_vector.push_back(k);

								}

							}

						}

					}
					sort_list.insert(pair<int, int>(counter, i));

				}

				//cout<<" Whats the fuck 1 "<<endl;
				if (SEQ_flag[i] != 1) {
					//fprintf(dereplication_name_file,"%d  ",i);
					for (int ik = 1; ik < strlen(SEQ_tag[i]); ik++)
						fprintf(dereplication_name_file, "%c", SEQ_tag[i][ik]);
					if (!temp_vector.empty())
						fprintf(dereplication_name_file, "__%d",
								((int) temp_vector.size() + 1));
					fprintf(dereplication_name_file, "   ");
					for (int ik = 1; ik < strlen(SEQ_tag[i]); ik++)
						fprintf(dereplication_name_file, "%c", SEQ_tag[i][ik]);
					if (!temp_vector.empty()) {
						fprintf(dereplication_name_file, ",");
						for (int tp = 0; tp < temp_vector.size(); tp++) {
							if (tp < temp_vector.size() - 1) {

								for (int ik = 1;
										ik < strlen(SEQ_tag[temp_vector[tp]]);
										ik++)
									fprintf(dereplication_name_file, "%c",
											SEQ_tag[temp_vector[tp]][ik]);
								fprintf(dereplication_name_file, ",");
							} else {
								for (int ik = 1;
										ik < strlen(SEQ_tag[temp_vector[tp]]);
										ik++)
									fprintf(dereplication_name_file, "%c",
											SEQ_tag[temp_vector[tp]][ik]);

								//fprintf(dereplication_name_file,"%d",temp_vector[tp]);
							}
						}
						//fprintf(dereplication_name_file,"\n");
						temp_vector.clear();
					}
					fprintf(dereplication_name_file, "\n");
				}

				//cout<<" Whats the fuck 2 "<<endl;
				//==========================================================
				//frintf();
				//==========================================================

				//count[i] = counter;
				//SEQ_flag[i] == 2;

			}
			fclose(dereplication_name_file);
			//===========================================================================

			cout << "Total number of sequences after Strict Dereplication "
					<< sort_list.size() << endl;

			int aA, cC, gG, tT;

			int ik = 0;

			int new_total = sort_list.size();

			//for(int ip = 0; ip < new_total; ip++)
			//{
			if ((frequency_list = new int[new_total]) == NULL) {
				printf("Memory Problem \n");
				help();
			}
			//}

			for (rmi = sort_list.rbegin(); rmi != sort_list.rend(); rmi++) {

				//cout << (*rmi).first << " " << (*rmi).second << endl;

				frequency_list[ik] = (*rmi).first;

				int len = SEQ_len[(*rmi).second];
				char *tmp_seq = SEQ_final[(*rmi).second];

				//cout<<tmp_seq[0]<<" * * * "<<tmp_seq[1]<<endl;

				int seq_tag_len = strlen(SEQ_tag[(*rmi).second]);

				char *tmp_seq_tag = SEQ_tag[(*rmi).second];

				for (int tl = 0; tl < seq_tag_len; tl++) {
					cd.push_back(tmp_seq_tag[tl]);
				}

				aA = 0;
				cC = 0;
				gG = 0;
				tT = 0;

				for (int kp = 0; kp < len; kp++) {

					//	cout<<tmp_seq[kp];
					//char ch = toupper(tmp[kp]);

					if (tmp_seq[kp] == 'A')

					{
						//count++;
						cd.push_back(tmp_seq[kp]);

						aA++;

					} else if (tmp_seq[kp] == 'C') {
						//count++;
						cd.push_back(tmp_seq[kp]);

						cC++;

					} else if (tmp_seq[kp] == 'G') {
						//count++;
						cd.push_back(tmp_seq[kp]);

						gG++;

					} else if (tmp_seq[kp] == 'T') {
						//count++;
						cd.push_back(tmp_seq[kp]);

						tT++;

					} else if (tmp_seq[kp] == 'N') {
						//count++;
						cd.push_back(tmp_seq[kp]);

					}

				}
				//cout<<endl;

				tmp.seq_id_size = seq_tag_len;
				tmpp.seq_id_size = seq_tag_len;

				globala[ik].push_back(aA);
				globala[ik].push_back(cC);
				globala[ik].push_back(gG);
				globala[ik].push_back(tT);

				global[ik] = cd;

				tmp.dy = 1;
				tmp.length = len;

				tmpp.dy = 1;

				dipsub[len][ik] = tmpp;

				dip[ik] = tmp;

				ik++;
				cd.clear();

			}

			//===========================Copy Part======================================

			//===========================================================================

			//================================ Release Memory ==================================

			sort_list.clear();
			delete[] SEQ_len;
			delete[] SEQ_flag;

			for (int i = 0; i < Total_sequence; i++) {
				delete[] SEQ_final[i];
				delete[] SEQ_tag[i];

			}
			//=================================================================

			printf("\nStrict Dereplication Ended \n");

			printf("\n Clustering started .... \n\n");

			if (kflag == 1)
				for (int i = kmin; i <= kmax;) {
					cout << " ---------------------------------- " << endl;
					cout << " \n Clustering at " << i << " Difference \n "
							<< endl;
					nbdf = i;
					cluster();
					if (i == kmax)
						break;
					i += ksteps;
					if (i > kmax)
						i = kmax;
				}
			else
				cluster();

			printf("\n Clustering ended .... \n");
			printf("\nElapsed time: %ld secs.\n", clock() / CLOCKS_PER_SEC);

			exit(10);

		}

		else {

			//nbdf = atoi(argv[1]);
			//strcpy(output, argv[3]);
			//strcpy(option, argv[4]);

			//===============================================================================================
			/*	strcpy(endgapflagg, argv[5]);
			 if (strcmp(endgapflagg, "--noendgaps") == 0)
			 endgapflag = 0;
			 else if (strcmp(endgapflagg, "--endgaps") == 0)
			 endgapflag = 1;
			 else {
			 printf("Please enter the correct option for endgaps \n");
			 printf(" --noendgaps    :For not counting the end gaps  \n");
			 printf(" --endgaps      :For counting the end gaps   \n");
			 exit(56);

			 }*/
			//=========================================================================================
			int aA, cC, gG, tT;

			map<short, short> a1;

			pFile = fopen(input_sequencefile, "r");

			if (pFile == NULL) {
				cout << " Input file not found .........   " << endl;
				help();
				//exit(24);
			}

			int count = 0;
			while (!feof(pFile)) {

				//int count = 0;
				do {
					fscanf(pFile, "%c", &ch);
					if ((ch != '>') && (ik == 0) && (count == 0)) {
						cout
								<< " Your input sequence fasta file is either empty or not formatted properly . . "
								<< endl;
						cout
								<< " Please check your fasta file whether it is empty or has some blank line in it . . "
								<< endl;
						cout << " Exiting the program . . ." << endl;
						exit(1);
					}
					if ((ch != 10) && (!isspace(ch))) {
						count++;
						cd.push_back(ch);

					}
				} while (ch != 10);

				if (!cd.empty()) {
					ik++;

					tmp.seq_id_size = count;
					tmpp.seq_id_size = count;

				}

				count = 0;

				aA = 0;
				cC = 0;
				gG = 0;
				tT = 0;

				do {
					fscanf(pFile, "%c", &ch);
					ch = toupper(ch);

					if (ch == 'A') {
						count++;
						cd.push_back(ch);

						aA++;

					} else if (ch == 'C') {
						count++;
						cd.push_back(ch);

						cC++;

					} else if (ch == 'G') {
						count++;
						cd.push_back(ch);

						gG++;

					} else if (ch == 'T') {
						count++;
						cd.push_back(ch);

						tT++;

					}

				} while ((!feof(pFile)) && (ch != '>'));

				if (!cd.empty()) {
					//ik++;
					globala[ik].push_back(aA);
					globala[ik].push_back(cC);
					globala[ik].push_back(gG);
					globala[ik].push_back(tT);

					global[ik] = cd;

					tmp.dy = 1;
					tmp.length = count;

					tmpp.dy = 1;

					dipsub[count][ik] = tmpp;

					dip[ik] = tmp;
				}
				cd.clear();
				if (ch == '>') {
					cd.push_back(ch);
					count = 1;
					//break;
				}

			}

			fclose(pFile); // ferme le flux et

			//===============================================

			//cout<<option<<endl;
			//===============================================

			if (kflag == 1)
				for (int i = kmin; i <= kmax;) {
					cout << " Clustering at " << i << " Difference " << endl;
					nbdf = i;
					cluster();
					if (i == kmax)
						break;
					i += ksteps;
					if (i > kmax)
						i = kmax;

					//i += int(kmax / ksteps);
				}
			else
				cluster();
		}

		printf("\nClustering ended .... \n\n");
		printf("\nElapsed time: %ld secs.\n", clock() / CLOCKS_PER_SEC);

	}

	//}

	return (0);

}

void cluster() {
	FILE* mm;
	vector<char> seed;
	vector<char> test;

	vector<char> seed_id;
	vector<char> test_id;

	vector<short> tseed;
	vector<short> ttest;

	map<short, vector<short> > ttseed;
	map<short, vector<short> > tttest;

	int length_seed;

	int range = 5;
	int counter = 0;
	map<int, seqsub> test_sub;
	map<int, seqsub>::iterator itsub;
	map<int, seqdet>::iterator it;
	//=========================================================================================

	map<int, map<int, seqsub> >::iterator ittdipsub;

	map<int, map<int, seqsub> >::reverse_iterator ritdipsub;

	ittdipsub = dipsub.begin();
	int shortest = (*ittdipsub).first;
	ritdipsub = dipsub.rbegin();
	int longest = (*ritdipsub).first;

	cout << " Length of shortest sequence  " << shortest << endl;
	cout << " Length of longest sequence   " << longest << endl;

	int middlerange = longest - shortest;

	//cout << " middle  " << middlerange << endl;

	//================================================================================
	map<int, map<int, seqsub> >::iterator itdipsub;
	//================================================================================

	int cluster = 0;
	int limit;
	if (strcmp(option, "--d_all") == 0) {
		rangeed = nbdf;
		if (nbdf <= 5) {
			limit = nbdf + nbdf + 1;
		} else
			limit = nbdf + 6;
	} else {

		limit = nbdf + 12;
		rangeed = nbdf + 15;
		//if (rangeed > middlerange)
		//rangeed = middlerange;

	}
	//==========================================================================================

	int sublooplower;
	int subloopupper;

	//==========================================================================================
	char buffer[50];
	char output_tmp[200];
	strcpy(output_tmp, output);
	//cout<<" ----------------- "<<output_tmp<<endl;
	sprintf(buffer, "_%d", nbdf);
	strcat(output_tmp, buffer);

	//==========================================================================================

	mm = fopen(output_tmp, "w");
	if (mm == NULL) {
		cout << "Error in opening outputfilef.fasta file " << endl;
		help();
		exit(1);
	}

	for (it = dip.begin(); it != dip.end(); it++) {

		if (argflag == 6)
			(dipsub[(it->second).length][it->first]).dy = 0;

		if ((it->second).dy != 0)

		{

			tseed = (globala[it->first]);

			seed = global[(it->first)];

			cluster++;

			fprintf(mm, "Cluster_%d\n", cluster);

			int seed_id_size = (it->second).seq_id_size;

			for (int i = 0; i < seed_id_size; i++) {
				fprintf(mm, "%c", seed[i]);
			}

			if (strict == 1) {
				fprintf(mm, "__%d", frequency_list[it->first]);
			}
			fprintf(mm, "\n");

			int seed_id_length = seed_id_size + (it->second).length;

			for (int i = seed_id_size; i < seed_id_length; i++) {

				fprintf(mm, "%c", seed[i]);
			}

			fprintf(mm, "\n");

			counter++;

			length_seed = (it->second).length;

			//================================================================================
			if (endgapflag == 0) {
				//cout<<"Hello "<<endl;
				sublooplower = abs(length_seed - shortest);
				subloopupper = abs(longest - length_seed);

			} else {
				sublooplower = rangeed;
				subloopupper = rangeed;
			}

			//========================================================================

			//cout<< sublooplower<<" "<<endgapflag<<" "<<limit<<endl;
			//cout<< subloopupper<<" "<<endgapflag<<" "<<limit<<endl;

			//================================================================================

			for (int iseek = length_seed - sublooplower;
					iseek <= length_seed + subloopupper; iseek++) {
				//=============================================================================================

				if (iseek < shortest)
					continue;
				else if (iseek > longest)
					break;
				//=============================================================================================

				test_sub = dipsub[iseek];

				if (!test_sub.empty()) {
					for (itsub = test_sub.begin(); itsub != test_sub.end();
							itsub++) {

						if ((itsub->second).dy != 0) {
							//================================================================================

							if (argflag == 6)
								ttest = globala[itsub->first];
							else
								ttest = globalb[itsub->first];

							//===================================================================
							if (endgapflag == 0) {

								limit = nbdf + 12 + abs(length_seed - iseek);
							}

							//================================================================================

							short distance_counter = 0;

							for (int ik = 0; ik < 4; ik++) {
								distance_counter = distance_counter
										+ abs(tseed[ik] - ttest[ik]);

								//cout<<"hello ---- > "<<distance_counter<<endl;

								if (distance_counter > limit)
									break;
							}

							if (distance_counter > limit)
								continue;

							else

							{
								//cout<<"hello ---- > "<<endl;

								int length_test = iseek;
								//================================================================================
								if (argflag == 6)
									test = global[(itsub->first)];
								else
									test = globall[(itsub->first)];
								//================================================================================
								int test_id_size = (itsub->second).seq_id_size;

								int test_id_length = length_test + test_id_size;
								//====================================================================================
								if (strcmp(option, "--d_all") == 0) {

									if (Distance(seed, test, seed_id_length,
											test_id_length, seed_id_size,
											test_id_size, range, nbdf)) {
										counter++;

										for (int i = 0; i < test_id_size; i++) {

											fprintf(mm, "%c", test[i]);

										}
										if (strict == 1) {
											fprintf(
													mm,
													"__%d",
													frequency_list[itsub->first]);
										}
										fprintf(mm, "\n");

										for (int i = test_id_size;
												i < test_id_length; i++) {

											fprintf(mm, "%c", test[i]);

										}

										fprintf(mm, "\n");

										(dipsub[length_test][itsub->first]).dy =
												0;
										//===================================================
										if (argflag == 6)
											(dip[itsub->first]).dy = 0;
										//===================================================

									}
								}

								//=====================================================================================
								else if (strcmp(option, "--d_pm") == 0) {

									if (Distancepm(seed, test, seed_id_length,
											test_id_length, seed_id_size,
											test_id_size, range, nbdf)) {
										counter++;

										for (int i = 0; i < test_id_size; i++) {

											fprintf(mm, "%c", test[i]);

										}

										fprintf(mm, "\n");

										for (int i = test_id_size;
												i < test_id_length; i++) {

											fprintf(mm, "%c", test[i]);

										}

										fprintf(mm, "\n");
										(dipsub[length_test][itsub->first]).dy =
												0;
										//===================================================
										if (argflag == 6)
											(dip[itsub->first]).dy = 0;
										//===================================================

									}
								}
								//=====================================================================================
								else if (strcmp(option, "--d_hl") == 0) {

									if (Distancehl(seed, test, seed_id_length,
											test_id_length, seed_id_size,
											test_id_size, range, nbdf)) {
										counter++;

										for (int i = 0; i < test_id_size; i++) {

											fprintf(mm, "%c", test[i]);

										}
										if (strict == 1) {
											fprintf(
													mm,
													"__%d",
													frequency_list[itsub->first]);
										}
										fprintf(mm, "\n");

										for (int i = test_id_size;
												i < test_id_length; i++) {

											fprintf(mm, "%c", test[i]);

										}

										fprintf(mm, "\n");
										(dipsub[length_test][itsub->first]).dy =
												0;
										//===================================================
										if (argflag == 6)
											(dip[itsub->first]).dy = 0;
										//===================================================

									}
								} else if (strcmp(option, "--d_hh") == 0) {

									if (Distancehh(seed, test, seed_id_length,
											test_id_length, seed_id_size,
											test_id_size, range, nbdf)) {
										counter++;

										for (int i = 0; i < test_id_size; i++) {

											fprintf(mm, "%c", test[i]);

										}

										fprintf(mm, "\n");

										for (int i = test_id_size;
												i < test_id_length; i++) {

											fprintf(mm, "%c", test[i]);

										}

										fprintf(mm, "\n");
										(dipsub[length_test][itsub->first]).dy =
												0;
										//===================================================
										if (argflag == 6)
											(dip[itsub->first]).dy = 0;
										//===================================================

									}
								}

								else {

									cout << "Please enter correct option..... "
											<< endl;
									cout << endl;
									help();

								}
								//=====================================================================================

							}

						}

					}
				}
			}

			fprintf(mm, "#_%d\n", counter);
			counter = 0;

		}

	}

	//===========================================================================

	if (argflag == 7) {
		fprintf(mm, "Unclustered_Sequences\n");
		int counter = 0;

		for (itdipsub = dipsub.begin(); itdipsub != dipsub.end(); itdipsub++) {
			test_sub = dipsub[itdipsub->first];
			int length_test = itdipsub->first;

			for (itsub = test_sub.begin(); itsub != test_sub.end(); itsub++) {

				if ((itsub->second).dy != 0) {

					test = globall[(itsub->first)];

					int test_id_size = (itsub->second).seq_id_size;

					int test_id_length = length_test + test_id_size;
					for (int i = 0; i < test_id_size; i++) {

						fprintf(mm, "%c", test[i]);

					}

					fprintf(mm, "\n");

					for (int i = test_id_size; i < test_id_length; i++) {

						fprintf(mm, "%c", test[i]);

					}
					fprintf(mm, "\n");
					counter++;

				}

				//*****************************Initialize for loop ******************************
				//(itsub->second).dy = 1;
				//***********************************************************

			}

		}
		fprintf(mm, "TotalUnclusteredSequence#_%d\n", counter);
	}
//********************************* Initialize for loop ***********************************
	for (it = dip.begin(); it != dip.end(); it++) {

		(it->second).dy = 1;
	}

	for (itdipsub = dipsub.begin(); itdipsub != dipsub.end(); itdipsub++) {
		for (itsub = dipsub[itdipsub->first].begin();
				itsub != dipsub[itdipsub->first].end(); itsub++)
			(itsub->second).dy = 1;
	}

	//********************************************************************

	//===========================================================================
	//if (argflag == 6)

	fclose(mm);
}

//================================================================================

//LOOP BEGINS HERE

//================================================================================

void clusterloop() {
	FILE* mm;
	vector<char> seed;
	vector<char> test;

	vector<char> seed_id;
	vector<char> test_id;

	vector<short> tseed;
	vector<short> ttest;

	map<short, vector<short> > ttseed;
	map<short, vector<short> > tttest;

	map<int, map<int, list<int> > > clustertrace;

	map<int, map<int, list<int> > >::iterator cltr;

	int length_seed;

	int range = 5;
	int counter = 0;
	map<int, seqsub> test_sub;

	map<int, seqdet>::iterator it;
	map<int, seqsub>::iterator itsub;

	//=======================================================================================

	map<int, map<int, seqsub> >::iterator ittdipsub;

	map<int, map<int, seqsub> >::reverse_iterator ritdipsub;

	ittdipsub = dipsub.begin();
	int shortest = (*ittdipsub).first;
	ritdipsub = dipsub.rbegin();
	int longest = (*ritdipsub).first;

	cout << " Length of shortest Sequence  " << shortest << endl;
	cout << " Length of longest  Sequence  " << longest << endl;

	int middlerange = (longest - shortest) / 2;

	//cout << " middle  " << middlerange << endl;

	//================================================================================
	map<int, map<int, seqsub> >::iterator itdipsub;
	//================================================================================

	int seqcounter = 0;
	//============================================================================
	int cluster = 0;
	int limit;
	//===================================================================

	limit = nbdf + 12;
	rangeed = nbdf + 15;

	//==========================================================================================

	int sublooplower;
	int subloopupper;

	//==========================================================================================

	//if (rangeed > middlerange)
	//rangeed = middlerange;

	//==================================================================

	//=======================================================================

	while ((seqcounter < seqnumber) && (nbdf < loop_upper_limit)) {

		//==========================================================================================

		cout << "Running at " << nbdf << " Difference " << endl;

		//list<int>seqlist;

		map<int, list<int> > looptrace;

		map<int, list<int> >::iterator lp;

		//==========================================================================================

		//==============================================
		for (it = dip.begin(); it != dip.end(); it++) {

			tseed = (globala[it->first]);

			seed = global[(it->first)];

			cluster++;

			int seed_id_size = (it->second).seq_id_size;

			int seed_id_length = seed_id_size + (it->second).length;

			counter++;

			length_seed = (it->second).length;
			//===================================================================================

			if (endgapflag == 0) {
				sublooplower = abs(length_seed - shortest);
				subloopupper = abs(longest - length_seed);

			} else {
				sublooplower = rangeed;
				subloopupper = rangeed;
			}

			//===============================================================================

			//cout<< sublooplower<<endl;
			//cout<< subloopupper<<endl;

			//=================================================================================

			for (int iseek = length_seed - sublooplower;
					iseek <= length_seed + subloopupper; iseek++) {

				//==========================================================================

				//cout<< iseek<<endl;

				if (iseek < shortest)
					continue;
				else if (iseek > longest)
					break;

				//==========================================================================

				test_sub = dipsub[iseek];
				if (!test_sub.empty()) {

					for (itsub = test_sub.begin(); itsub != test_sub.end();
							itsub++) {

						if ((itsub->second).dy != 0) {
							//================================================================================

							if (endgapflag == 0)
								limit = nbdf + 12 + abs(length_seed - iseek);

							ttest = globalb[itsub->first];

							//================================================================================

							short distance_counter = 0;

							for (int ik = 0; ik < 4; ik++) {
								distance_counter = distance_counter
										+ abs(tseed[ik] - ttest[ik]);

								if (distance_counter > limit)
									break;
							}

							if (distance_counter > limit)
								continue;

							else

							{
								//cout<<"hello "<<endl;

								int length_test = iseek;
								//================================================================================

								test = globall[(itsub->first)];
								//================================================================================
								int test_id_size = (itsub->second).seq_id_size;

								int test_id_length = length_test + test_id_size;
								//====================================================================================

								//=====================================================================================

								if (Distancehl(seed, test, seed_id_length,
										test_id_length, seed_id_size,
										test_id_size, range, nbdf)) {

									//cout<<"hello 0000 "<<endl;

									clustertrace[cluster][iseek].push_back(
											itsub->first);

									(dipsub[length_test][itsub->first]).distance =
											nbdf;

									//==========================================================

									looptrace[iseek].push_back(itsub->first);

									//==========================================================
									//==========================================================

								}

								//=====================================================================================

							}

						}

					}

				} else
					continue;

				//===============================================================================

				//===================================================================================

			}

		}

		//============================================================================

		if (looptrace.size() > 0) {

			//cout<<"helloooo "<<endl;

			for (lp = looptrace.begin(); lp != looptrace.end(); lp++) {
				list<int> tmpp;
				list<int>::iterator tpp;
				tmpp = lp->second;
				tmpp.sort();
				tmpp.unique();

				int size = tmpp.size();
				seqcounter = seqcounter + tmpp.size();

				for (tpp = tmpp.begin(); tpp != tmpp.end(); tpp++) {

					(dipsub[lp->first][*tpp]).dy = 0;

				}

			}
		}

		//=============================================================================

		nbdf++;
		limit = nbdf + 12;
		rangeed = nbdf + 15;
		cluster = 0;

	}

	//===========================================================================

	mm = fopen(output, "w");
	if (mm == NULL) {
		cout << "Error in opening output file " << endl;
		exit(1);
	}

	for (it = dip.begin(); it != dip.end(); it++) {
		cluster++;
		counter = 0;

		seed = global[(it->first)];

		int seed_id_size = (it->second).seq_id_size;

		int seed_id_length = seed_id_size + (it->second).length;

		fprintf(mm, "Cluster_%d\n", cluster);

		for (int i = 0; i < seed_id_size; i++) {
			fprintf(mm, "%c", seed[i]);
		}
		fprintf(mm, "\n");

		for (int i = seed_id_size; i < seed_id_length; i++) {

			fprintf(mm, "%c", seed[i]);
		}

		fprintf(mm, "\n");

		//============================================================================

		counter++;

		//============================================================================

		length_seed = (it->second).length;

		//============================================================================

		//============================================================================

		for (int iseek = shortest; iseek <= longest; iseek++) {

			//============================================================================

			list<int> tempp = clustertrace[cluster][iseek];

			tempp.sort();

			list<int>::iterator temppp;

			test_sub = dipsub[iseek];

			int length_test = iseek;

			if (!test_sub.empty()) {

				for (temppp = tempp.begin(); temppp != tempp.end(); temppp++) {

					int test_id_size = test_sub[*temppp].seq_id_size;

					int test_id_length = length_test + test_id_size;

					counter++;
					test = globall[*temppp];
					for (int i = 0; i < test_id_size; i++) {

						fprintf(mm, "%c", test[i]);

					}
					fprintf(mm, "##Distance_<=_%d", test_sub[*temppp].distance);

					fprintf(mm, "\n");

					for (int i = test_id_size; i < test_id_length; i++) {

						fprintf(mm, "%c", test[i]);

					}

					fprintf(mm, "\n");

				}

			}

		}

		fprintf(mm, "#_%d\n", counter);

	}

	//===========================================================================

	if (seqcounter < seqnumber) {

		//==========================================================================

		//==========================================================================

		fprintf(mm, "Unclustered_Sequences\n");
		int counter = 0;

		for (itdipsub = dipsub.begin(); itdipsub != dipsub.end(); itdipsub++) {
			test_sub = dipsub[itdipsub->first];
			int length_test = itdipsub->first;

			for (itsub = test_sub.begin(); itsub != test_sub.end(); itsub++) {

				if ((itsub->second).dy != 0) {

					test = globall[(itsub->first)];

					int test_id_size = (itsub->second).seq_id_size;

					int test_id_length = length_test + test_id_size;
					for (int i = 0; i < test_id_size; i++) {

						fprintf(mm, "%c", test[i]);

					}

					fprintf(mm, "\n");

					for (int i = test_id_size; i < test_id_length; i++) {

						fprintf(mm, "%c", test[i]);

					}
					fprintf(mm, "\n");
					counter++;

				}

			}

		}
		fprintf(mm, "TotalUnclusteredSequence#_%d\n", counter);
	}
	//===========================================================================
	//if (argflag == 6)

	fclose(mm);
}

//================================================================================

bool Distance(vector<char>& s1, vector<char>& s2, int s1_length, int s2_length,
		int seed_id_size, int test_id_size, int msearch, int nbdf)

		{

	int k = (s1_length - seed_id_size) - (s2_length - test_id_size);

	float dist = 0;

	int c = 0;

	int ora = seed_id_size;

	int orb = test_id_size;

	int common = 0;

	if (nbdf == 1) {
		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{

			//cout << "hi  " << endl;

			//cout << s1[c + ora] << " " << s2[c + orb] << " " << k << nbdf
			//	<< endl;

			if (endgapflag == 0) {

				if ((c + ora == s1_length - 1) || (c + orb == s2_length - 1))
					return (1);

			}

			else if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
				return (1);

			if (s1[c + ora] == s2[c + orb]) {

				common++;
				c++;
			} else

			{
				dist++;
				if (dist > nbdf) {
					return (0);
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((k == 0) && (c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						ora = ora + i;
						orb = orb + i;
						if (i > 1) {
							dist = dist + i - 1;
							if (dist > nbdf) {
								return (0);
							}
						}

						common++;

						break;

					}

					else if ((k < 0) && (c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						orb = orb + i;

						common++;
						if (i > 1) {
							dist = dist + i - 1;
							if (dist > nbdf) {
								return (0);
							}
						}

						break;

					}

					else if ((k > 0) && (c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{

						ora = ora + i;

						if (i > 1) {
							dist = dist + i - 1;
							if (dist > nbdf) {
								return (0);
							}
						}
						common++;

						break;

					}

				}
				c++;

			}

		}

		if (endgapflag == 1) {

			if ((ceil(
					(double) (s1_length - seed_id_size + s2_length
							- test_id_size) / 2) - common) <= nbdf) {

				return (1);
			} else {

				return (0);
			}
		} else if (endgapflag == 0) {

			if ((c + ora == s1_length) || (c + orb == s2_length))
				return (1);
		} else {

			return (0);
		}

	}

	else if (nbdf == 0) {
		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{

			//cout<<"hi  "<<endl;
			if (s1[c + ora] != s2[c + orb]) {

				return (0);

			} else

				c++;
		}

		if (endgapflag == 1) {

			if ((ceil(
					(double) (s1_length - seed_id_size + s2_length
							- test_id_size) / 2) - c) <= nbdf) {

				return (1);
			} else if (((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
					|| ((c + ora == s1_length) && (c + orb == s2_length)))
				return (1);
			else
				return (0);
		} else if (endgapflag == 0) {

			if ((c + ora == s1_length - 1) || (c + orb == s2_length - 1)
					|| (c + ora == s1_length) || (c + orb == s2_length))
				return (1);
			else
				return (0);
		}

		//return (1);
		//========================================================================================

		// No end gap flag

		//========================================================================================

	}

	else {

		if (Distancemain(s1, s2, s1_length, s2_length, msearch, nbdf, 0,
				seed_id_size, test_id_size))
			return (1);
		else
			return (0);
	}

}

bool Distancemain(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb) {

	signed int kk = (s1_length - c - ora) - (s2_length - c - orb);

	if (endgapflag == 1) {
		if (((abs)(kk) <= nbdf)
				&& ((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {
			return (1);
		}
	} else if (endgapflag == 0) {
		if ((c + ora == s1_length - 1) || (c + orb == s2_length - 1)) {
			return (1);
		}
	}

	if (nbdf <= 1) {

		if (Distanceone(s1, s2, s1_length, s2_length, msearch, nbdf, c, ora,
				orb))
			return (1);
		else
			return (0);
	}

	int dist = 0;

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{

		//================================================================
		//cout << s1[c + ora] << " " << s2[c + orb] << " " << kk << nbdf << endl;
		//================================================================
		if (endgapflag == 1) {
			if (((abs)(kk) <= nbdf)
					&& ((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {
				return (1);
			}
		} else if (endgapflag == 0) {
			if ((c + ora == s1_length - 1) || (c + orb == s2_length - 1)) {
				return (1);
			}
		}

		if (s1[c + ora] == s2[c + orb]) {

			if (endgapflag == 1) {
				if (((abs)(kk) <= nbdf)
						&& ((c + ora == s1_length - 1)
								|| (c + orb == s2_length - 1))) {
					return (1);
				}
			} else if (endgapflag == 0) {
				if ((c + ora == s1_length - 1) || (c + orb == s2_length - 1)) {
					return (1);
				}

			}
			c++;

		}

		else

		{

			dist++;
			signed int k = (s1_length - (c + ora)) - (s2_length - (c + orb));

			if (k < 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int orbb = orb + i;

						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);

						//===================================================

						//cout<<"Hellooo <<< "<<distt<<" "<<endgapflag<<endl;

						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;
						//}

						//====================================================

						else {

							//cout<<"Hellooo 1 ==== "<<distt<<endl;

							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);

							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{

						int oraa = ora + i;
						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);

						//=================================================
						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;
						//================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						int oraa = ora + i;
						int orbb = orb + i;
						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						//====================================================
						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;
						//====================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}
				return (0);

			}

			else if (k == 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						int oraa = ora + i;
						int orbb = orb + i;
						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);

						//===================================================
						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;

						//===================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;

						}

					}

				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{

						int oraa = ora + i;
						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						//===================================================
						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;
						//====================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;

						}

					}
				}

				for (int i = 1; i <= msearch; i++) {

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int orbb = orb + i;

						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						//===================================================
						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;
						//=====================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);
							} else
								break;

						}

					}
				}
				return (0);

			}

			else if (k > 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{

						int oraa = ora + i;
						int distt = dist + i - 1;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);

						//====================================================

						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;
						//===================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						int oraa = ora + i;
						int orbb = orb + i;
						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);

						//================================================
						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 1) && ((distt > nbdf)))
							break;
						//=================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}

				for (int i = 1; i <= msearch; i++) {

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int orbb = orb + i;

						int distt = dist + i - 1;
						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);

						//=================================================

						//if (endgapflag == 1) {
						if ((endgapflag == 1)
								&& (((abs)(kk) > nbdf - distt) || (distt > nbdf)))
							break;
						//} else if (endgapflag == 0) {
						else if ((endgapflag == 0) && ((distt > nbdf)))
							break;

						//==================================================

						else {
							if (Distancemain(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);
							} else
								break;
						}

					}
				}
				return (0);

			}

		}

		if (endgapflag == 1) {
			if (((abs)(kk) <= nbdf)
					&& ((c + ora == s1_length) || (c + orb == s2_length))) {
				return (1);
			}
		} else if (endgapflag == 0) {
			if ((c + ora == s1_length) || (c + orb == s2_length)) {
				return (1);
			}
		}

	}

}

//==================================================================================================
bool Distanceone(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb)

		{

	int s11_length = s1_length - c - ora;
	int s22_length = s2_length - c - orb;

	int k = (s11_length - s22_length);

	float dist = 0;

	int common = 0;

	//cout<<"I am here "<<endl;

	if (endgapflag == 1) {
		if (((abs)(k) <= nbdf)
				&& (((c + ora == s1_length) && (c + orb == s2_length))
						|| ((c + ora == s1_length - 1)
								&& (c + orb == s2_length - 1)))) {
			return (1);
		}
	} else if (endgapflag == 0) {
		if ((c + ora == s1_length - 1) || (c + orb == s2_length - 1)
				|| (c + ora == s1_length) || (c + orb == s2_length)) {
			return (1);
		}
	}

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{
		if (endgapflag == 1) {
			if ((c == s1_length - ora - 1) && (c == s2_length - orb - 1)) {
				return (1);
			}
		} else if (endgapflag == 0) {
			if ((c == s1_length - ora - 1) || (c == s2_length - orb - 1)) {
				return (1);
			}

		}

		if (s1[c + ora] == s2[c + orb]) {

			if (endgapflag == 1) {
				if ((c == s1_length - ora - 1) && (c == s2_length - orb - 1)) {
					return (1);
				}
			} else if (endgapflag == 0) {
				if ((c == s1_length - ora - 1) || (c == s2_length - orb - 1)) {
					return (1);
				}

			}

			common++;
			c++;

		} else

		{
			dist++;
			if (dist > nbdf) {
				return (0);
			}

			for (int i = 1; i <= msearch; i++)

			{

				if ((k == 0) && (c + ora + i < s1_length)

				&& (s1[c + ora + i] == s2[c + orb + i]))

				{

					ora = ora + i;
					orb = orb + i;
					common++;
					if (i > 1) {
						dist = dist + i - 1;
						if (dist > nbdf) {
							return (0);
						}
					}

					break;

				}

				else if ((k < 0) && (c + orb + i < s2_length)

				&& (s1[c + ora] == s2[c + orb + i]))

				{

					orb = orb + i;

					common++;
					if (i > 1) {
						dist = dist + i - 1;
						if (dist > nbdf) {
							return (0);
						}
					}

					break;

				}

				else if ((k > 0) && (c + ora + i < s1_length)
						&& (s1[c + ora + i] == s2[c + orb]))

						{

					ora = ora + i;
					common++;

					if (i > 1) {
						dist = dist + i - 1;
						if (dist > nbdf) {
							return (0);
						}
					}

					break;

				}

			}
			c++;

		}
		if (endgapflag == 1) {
			if ((c + ora == s1_length) && (c + orb == s2_length)) {
				return (1);
			}
		} else if (endgapflag == 0) {
			if ((c + ora == s1_length) || (c + orb == s2_length)) {
				return (1);
			}
		}

	}

	if (endgapflag == 1) {

		if ((ceil((double) (s11_length + s22_length) / 2) - common) <= nbdf) {

			return (1);
		} else {

			return (0);
		}
	} else if (endgapflag == 0) {

		if ((c + ora == s1_length) || (c + orb == s2_length))
			return (1);
	}

}
//=====================================================================================

//===================================== PointMutation ================================

bool Distancepm(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int seed_id_size, int test_id_size, int msearch,
		int nbdf)

		{

	int k = (s1_length - seed_id_size) - (s2_length - test_id_size);

	int c = 0;

	int ora = seed_id_size;

	int orb = test_id_size;

	int common = 0;

	float gdist = 0;

	if (nbdf == 1) {
		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{

			if (s1[c + ora] == s2[c + orb]) {

				common++;
				c++;
			} else

			{

				for (int i = 1; i <= msearch; i++)

				{

					if ((k == 0) && (c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						float dist = 0;
						ora = ora + i;
						orb = orb + i;

						if (i > 1) {
							dist = i - 1 + gdist;
						} else
							dist = i + gdist;
						gdist = dist;
						if (gdist > nbdf) {
							return (0);
						}

						common++;

						break;

					}

					else if ((k < 0) && (c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						orb = orb + i;

						common++;

						break;

					}

					else if ((k > 0) && (c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{

						ora = ora + i;

						common++;

						break;

					}
					if (i == msearch)
						gdist++;
				}
				c++;

			}

		}

		return (1);

	}
	//======================================================NBDF==0======================
	else if (nbdf == 0) {
		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{

			if (s1[c + ora] == s2[c + orb]) {

				common++;
				c++;
			} else

			{

				for (int i = 1; i <= msearch; i++)

				{

					if ((k == 0) && (c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						float dist = 0;
						ora = ora + i;
						orb = orb + i;

						if (i > 1) {
							dist = i - 1;
						} else
							dist = i;
						if (dist > nbdf) {
							return (0);
						}

						common++;

						break;

					}

					else if ((k < 0) && (c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						orb = orb + i;

						common++;

						break;

					}

					else if ((k > 0) && (c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{

						ora = ora + i;

						common++;

						break;

					}
					if (i == msearch)
						return (0);
				}

				c++;

			}

		}

		return (1);

	}

	else {

		if (Distancemainpm(s1, s2, s1_length, s2_length, msearch, nbdf, 0,
				seed_id_size, test_id_size))
			return (1);
		else
			return (0);
	}

}

bool Distancemainpm(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb) {

	if (nbdf <= 1) {

		if (Distanceonepm(s1, s2, s1_length, s2_length, msearch, nbdf, c, ora,
				orb))
			return (1);
		else
			return (0);
	}

	signed int kk = (s1_length - c - ora) - (s2_length - c - orb);

	if (((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {
		return (1);
	}

	int dist = 0;

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{

		if (s1[c + ora] == s2[c + orb]) {

			if (((abs)(kk) <= nbdf)
					&& ((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {
				return (1);
			}

			c++;

		}

		else

		{

			signed int k = (s1_length - (c + ora)) - (s2_length - (c + orb));

			if (k < 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int orbb = orb + i;

						int distt = dist;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);

							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{

						int oraa = ora + i;

						int distt = dist;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						int distt = 0;
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = i - 1;
						} else
							distt = i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}
				return (0);

			}

			else if (k == 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						int distt = 0;
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = i - 1;
						} else
							distt = i;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;

						}

					}

				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{

						int oraa = ora + i;

						int distt = dist;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;

						}

					}
				}

				for (int i = 1; i <= msearch; i++) {

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int orbb = orb + i;

						int distt = dist;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);
							} else
								break;

						}

					}
				}
				return (0);

			}

			else if (k > 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{

						int oraa = ora + i;

						int distt = dist;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						int distt = 0;
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = i - 1;
						} else
							distt = i;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}

				for (int i = 1; i <= msearch; i++) {

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int orbb = orb + i;

						int distt = dist;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainpm(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);
							} else
								break;
						}

					}
				}
				return (0);

			}

		}

	}

}

bool Distanceonepm(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb)

		{

	int s11_length = s1_length - c - ora;
	int s22_length = s2_length - c - orb;

	int k = (s11_length - s22_length);

	int common = 0;
	float gdist = 0;

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{

		if (s1[c + ora] == s2[c + orb]) {

			if ((c == s1_length - ora - 1) || (c == s2_length - orb - 1)) {

				return (1);
			}

			common++;
			c++;

		} else

		{
			for (int i = 1; i <= msearch; i++)

			{

				if ((k == 0) && (c + ora + i < s1_length)

				&& (s1[c + ora + i] == s2[c + orb + i]))

				{
					float dist = 0;
					ora = ora + i;
					orb = orb + i;
					common++;

					if (i > 1) {
						dist = i - 1 + gdist;
					} else
						dist = i + gdist;
					gdist = dist;
					if (gdist > nbdf) {

						return (0);

					}

					break;

				}

				else if ((k < 0) && (c + orb + i < s2_length)

				&& (s1[c + ora] == s2[c + orb + i]))

				{

					orb = orb + i;
					common++;
					break;

				}

				else if ((k > 0) && (c + ora + i < s1_length)
						&& (s1[c + ora + i] == s2[c + orb]))

						{

					ora = ora + i;
					common++;
					break;

				}
				if (i == msearch)
					gdist++;
			}
			c++;

		}

	}

	return (1);

}

//========================================== PointMutation+HomoPolymer+HL+Examples(AA, CC, TT)======================

//==================================================================================================================

bool Distancehl(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int seed_id_size, int test_id_size, int msearch,
		int nbdf)

		{

	int k = (s1_length - seed_id_size) - (s2_length - test_id_size);

	//======================================================================================

	/*
	 *
	 * if (endgapflag == 1)
	 {
	 if ( abs(k) > nbdf + 6)
	 return (0);

	 }
	 */
	//======================================================================================
	//cout<<"value of k "<<k<<endl;
	int c = 0;

	int ora = seed_id_size;

	int orb = test_id_size;

	int flagi = 0;

	int flagd = 0;

	int common = 0;

	int flag = 0;

	int counterr = 0;

	int gdist = 0;

	if (nbdf == 1) {

		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{
			//===================================================================
			//cout << "  " << s1[c + ora] << " " << s2[c + orb] << " "
			//	<< endgapflag << endl;
			//===================================================================

			if (c + orb >= s2_length - 1) {

				if (endgapflag == 1) {

					int endgap = s1_length - (c + ora);
					if (endgap > nbdf)
						return (0);
				}

			} else if (c + ora >= s1_length - 1) {

				if (endgapflag == 1) {

					int endgap = s2_length - (c + orb);
					if (endgap > nbdf)
						return (0);
				}

			}

			if (s1[c + ora] == s2[c + orb]) {
				if (c >= 1) {
					if (s1[c + ora] == s1[c + ora - 1]) {
						if (s1[c + ora - 1] == s2[c + orb - 1])
							counterr++;
						else
							counterr = 1;
					} else
						counterr = 1;
				} else
					counterr = 1;

				//================================================================================
				if (endgapflag == 0) {
					if ((c + ora == s1_length - 1)
							|| (c + orb == s2_length - 1))
						return (1);
				}
				//==========================================================================================
				if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
					return (1);

				else if ((c + ora < s1_length) && (c + orb == s2_length - 1)) {

					//=====================================================================================

					//====================================================================================

					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);
					else if (flagint == 1) {
						for (int ik = c + ora + 1; ik < s1_length; ik++)
							if (s1[ik] == s2[c + orb])
								flagint = 1;
							else {
								gdist++;
								flagint = 0;
								break;
							}
					}

					if ((flagint == 1) && (gdist <= 1))
						return (1);
					else
						return (0);

				} else if ((c + ora == s1_length - 1)
						&& (c + orb < s2_length)) {

					//=========================================================================

					//=========================================================================

					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);

					if (flagint == 1) {
						for (int ik = c + orb + 1; ik < s2_length; ik++)
							if (s1[c + ora] == s2[ik])
								flagint = 1;
							else {
								gdist++;
								flagint = 0;
								break;
							}
					}
					if ((flagint == 1) && (gdist <= 1))
						return (1);
					else
						return (0);
				}

				//==========================================================================================

				common++;
				c++;

			} else

			{
				//===============================================================
				//cout << "  **** " << s1[c + ora] << " " << s2[c + orb] << endl;
				//===============================================================
				if (counterr >= 2) {
					flag = 1;

				} else if (counterr == 1)
					flag = 2;
				counterr = 0;
				if ((flag == 1) && (s1[c + ora] == s1[c + ora - 1])) {
					flagi = 1;
					flagd = 0;
					flag = 0;
				} else if ((flag == 1) && (s2[c + orb] == s2[c + orb - 1])) {
					flagd = 1;
					flagi = 0;
					flag = 0;
				} else if ((flag == 2) && (s1[c + ora] == s1[c + ora - 1])) {
					flagi = 1;
					flagd = 0;
					flag = 0;
				} else if ((flag == 2) && (s2[c + orb] == s2[c + orb - 1])) {
					//if (s1[c + ora] == s1[c + ora - 1]) {

					flagd = 1;
					flagi = 0;
					flag = 0;
					/*} else {
					 flagi = 0;
					 flagd = 0;
					 flag = 0;
					 }*/
				} else {
					flagi = 0;
					flagd = 0;
					flag = 0;
				}
				for (int i = 1; i <= msearch; i++)

				{

					if ((k == 0) && (c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						//cout<<"Value of flagi "<<flagi<<i<<endl;

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist += 0;
							} else
								gdist += 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist += 0;
							} else
								gdist += 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist += dist;
						}

						//cout<<"Value of gdist "<<gdist<<i<<endl;

						counterr++;
						ora = ora + i;
						orb = orb + i;

						if (gdist > nbdf) {
							return (0);
							//break;
						}

						common++;

						break;

					}

					else if ((k < 0) && (c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						counterr++;
						float dist = 1;

						/*if (flagi == 1) {
						 dist = 0;
						 }

						 else {
						 if (i > 1) {
						 dist = dist + i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;
						 }
						 if (gdist > nbdf) {
						 return (0);
						 }*/
						//==========================================================================================
						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist++;
						} else if ((i == 1) && (flagd == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;

						}
						gdist = dist;
						if ((gdist > nbdf)
								|| ((s1[c + ora + 1] != s2[c + orb + i + 1]))) {
							//cout<<"Gdist __ *** "<<gdist<<endl;

							goto levelD;
							//return (0);
						}
						//====================================================================================

						orb = orb + i;

						common++;

						break;

					}

					else if ((k > 0) && (c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{
						counterr++;
						float dist = 1;
						//==========================================================================================
						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s1[c + ora] == s1[c + ora + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist = i - 1;

						} else if ((i == 1) && (flagi == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;

							} else {
								dist = i + gdist;

							}
						}
						gdist = dist;

						//cout<<" I Flagi "<<i<<" "<<flagi<<s1[c + ora + i]<<" "<<s2[c + orb]<<" "<<gdist<<endl;
						//cout<<"Gdist "<<gdist<<endl;

						if ((gdist > nbdf)
								|| (s1[c + ora + i + 1] != s2[c + orb + 1])) {
							//cout<<"Gdist __ "<<i<<s1[c + ora+i]<<" "<<s2[c + orb+i]<<" "<<gdist<<endl;
							goto levelD;
							//return (0);
						}

						//==========================================================================================

						ora = ora + i;

						common++;

						break;

					}

					//===========================================================
					levelD: if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						//cout << "Hello 1 : " << gdist << endl;

						//cout << s2[c + orb] << " " << s2[c + orb + i] << endl;

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist += 0;
							} else
								gdist += 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist += 0;
							} else
								gdist += 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist += dist;
						}

						//cout << "Hello 2: " << gdist << endl;

						counterr++;
						ora = ora + i;
						orb = orb + i;

						if (gdist > nbdf) {
							return (0);
							//break;
						}
						//cout << "Hello " << gdist << endl;
						common++;

						break;

					}

					//===========================================================
					if (i == msearch)
						gdist++;
					if (gdist > nbdf) {
						return (0);
					}

				}
				c++;

			}

		}

		if ((c + ora == s1_length) || (c + orb == s2_length))
			return (1);
		else
			return (0);

	}

	//====================================[ NBDF==0 ]================================================

	else if (nbdf == 0) {

		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{

			if (c + orb >= s2_length - 1) {

				if (endgapflag == 1) {

					int endgap = s1_length - (c + ora);
					if (endgap > nbdf + 1)
						return (0);
				}

			} else if (c + ora >= s1_length - 1) {

				if (endgapflag == 1) {

					int endgap = s2_length - (c + orb);
					if (endgap > nbdf + 1)
						return (0);
				}

			}

			if (s1[c + ora] == s2[c + orb]) {
				if (c >= 1) {
					if (s1[c + ora] == s1[c + ora - 1]) {
						if (s1[c + ora - 1] == s2[c + orb - 1])
							counterr++;
						else
							counterr = 1;

					} else
						counterr = 1;
				} else
					counterr = 1;
				//=================================================================================
				if (endgapflag == 0) {
					if ((c + ora == s1_length - 1)
							|| (c + orb == s2_length - 1))
						return (1);
				}

				//==============================================================================if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))

				//cout << "=====  " << s1[c + ora] << " " << s2[c + orb] << " "
				//	<< endgapflag << endl;

				//==========================================================================================
				if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
					return (1);

				else if ((c + ora < s1_length) && (c + orb == s2_length - 1)) {

					//=======================================================================================

					//=======================================================================================

					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);
					if (flagint == 1) {
						for (int ik = c + ora + 1; ik < s1_length; ik++)
							if (s1[ik] == s2[c + orb])
								flagint = 1;
							else
								flagint = 0;
					}

					if (flagint == 1)
						return (1);
					else
						return (0);

				} else if ((c + ora == s1_length - 1)
						&& (c + orb < s2_length)) {
					//=======================================================================================

					//=======================================================================================

					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);

					if (flagint == 1) {
						for (int ik = c + orb + 1; ik < s2_length; ik++)
							if (s1[c + ora] == s2[ik])
								flagint = 1;
							else
								flagint = 0;
					}
					if (flagint == 1)
						return (1);
					else if (flagint == 0)
						return (0);
				}

				//==========================================================================================
				common++;
				c++;

			} else

			{

				//============================================================

				//cout << " ___ " << s1[c + ora] << " " << c + ora << " "
				//	<< s1_length << " " << c + orb << " " << s2_length
				//<< " " << s2[c + orb] << " " << endgapflag << flag
				//<< endl;

				//if(c+ora == s1)

				//============================================================

				if (counterr >= 2) {
					flag = 1;

				} else if (counterr == 1)
					flag = 2;
				counterr = 0;
				if ((flag == 1) && (s1[c + ora] == s1[c + ora - 1])) {
					flagi = 1;
					flagd = 0;
					flag = 0;
				} else if ((flag == 1) && (s2[c + orb] == s2[c + orb - 1])) {
					flagd = 1;
					flagi = 0;
					flag = 0;
				} else if ((flag == 2) && (s1[c + ora] == s1[c + ora - 1])) {
					flagi = 1;
					flagd = 0;
					flag = 0;
				} else if ((flag == 2) && (s2[c + orb] == s2[c + orb - 1])) {
					//if (s1[c + ora] == s1[c + ora - 1]) {

					flagd = 1;
					flagi = 0;
					flag = 0;
					/*} else {
					 flagi = 0;
					 flagd = 0;
					 flag = 0;
					 } */
				} else {
					flagi = 0;
					flagd = 0;
					flag = 0;
				}

				for (int i = 1; i <= msearch; i++)

				{
					//cout << "Value of i " << i << s1[c + ora + i] << " "
					//	<< s2[c + orb + i] << endl;

					//cout << "Value of i " << c + ora + i << " " << c + orb + i
					//	<< " " << s1_length << k << endl;

					if ((k == 0) && (c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						//cout << "Value of i <<<>>> " << i << s1[c + ora + i]
						//	<< " " << s2[c + orb + i] << endl;

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist = 0;
							} else
								gdist = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist = 0;
							} else
								gdist = 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist = dist;
						}

						//cout << " value of gdist " << gdist << " " << i << " "
						//	<< flagi << endl;

						counterr++;
						//float dist = 1;
						ora = ora + i;
						orb = orb + i;
						/*if (i > 1) {
						 dist = dist + i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;*/

						if (gdist > nbdf) {
							return (0);
						}

						common++;

						break;

					}

					//else if
					//===========================================

					//===========================================

					else if ((k < 0) && (c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{
						//cout<<"value of i "<<i<<" "<<flagd<<" "<<s1[c+ora]<<" "<<s2[c+orb+i]<<endl;

						counterr++;
						float dist = 1;

						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist = i - 1;
						} else if ((i == 1) && (flagd == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;

						}
						gdist = dist;
						if ((gdist > nbdf)
								|| (s1[c + ora + 1] != s2[c + orb + i + 1])) {
							goto levelf;
							//return (0);
						}

						//====================================================================================

						orb = orb + i;

						common++;

						break;

					}

					//=========================================================
					//else if

					//=========================================================

					else if ((k > 0) && (c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{

						counterr++;
						float dist = 1;

						//	cout<<"value of i "<<i<<" "<<flagi<<" "<<s1[c+ora + i]<<" "<<s2[c+orb]<<endl;

						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								//cout<<"value of i i "<<i<<" "<<ik<<" "<<flagi<<" "<<s1[c+ora]<<" "<<s1[c+ora+ik]<<endl;
								if (s1[c + ora] == s1[c + ora + ik]) {
									flagdin = 1;
								} else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist = i - 1;
							//cout<<"dist  "<<dist<<endl;

						} else if ((i == 1) && (flagi == 1)) {
							//==========================================================================================
							//if(s1[c + ora] == s1[c + ora + i])
							dist = 0;
							//else
							//dist = i + gdist;
							//==========================================================================================
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;

							} else {
								dist = i + gdist;

							}
							//gdist = dist;
						}
						gdist = dist;
						if ((gdist > nbdf)
								|| (s1[c + ora + i + 1] != s2[c + orb + 1])) {
							goto levelf;
							//return (0);
						}

						ora = ora + i;

						common++;

						break;

					}
					levelf: if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						//cout << "Hello " << gdist << endl;

						//cout << s2[c + orb] << " " << s2[c + orb + i] << endl;

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist = 0;
							} else
								gdist = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist = 0;
							} else
								gdist = 1;
						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist = dist;
						}

						counterr++;
						ora = ora + i;
						orb = orb + i;

						if (gdist > nbdf) {
							return (0);
							//break;
						}

						common++;

						break;

					}

					if (i == msearch)
						return (0);
				}

				c++;

			}
			//============================================================================================
		}

		if ((c + ora == s1_length) || (c + orb == s2_length))
			return (1);
		else
			return (0);

	}

	else {
		globalflag = 1;
		if (Distancemainhl(s1, s2, s1_length, s2_length, msearch, nbdf, 0,
				seed_id_size, test_id_size))
			return (1);
		else
			return (0);
	}

}

bool Distancemainhl(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb) {

	int counterr;
	if (globalflag == 1) {
		counterr = 0;
		globalflag = 0;
	} else
		counterr = 1;

	int flagi = 0;
	int flagd = 0;

	//===================================================================================================

	if (endgapflag == 0) {
		if ((c + ora == s1_length) || (c + orb == s2_length))
			return (1);
	}

	else if (endgapflag == 1) {
		if ((c + ora == s1_length) && (c + orb == s2_length))
			return (1);
	}

	//=========================================================================================================
	if (nbdf <= 1) {

		if (Distanceonehl(s1, s2, s1_length, s2_length, msearch, nbdf, c, ora,
				orb))
			return (1);
		else
			return (0);
	}

	signed int kk = (s1_length - c - ora) - (s2_length - c - orb);

	//=========================================================================================
	if (endgapflag == 1) {

		if ((c + orb >= s2_length - 1)) {

			int endgap = s1_length - (c + ora);
			if (endgap > nbdf)
				return (0);

		}

		else if (c + ora >= s1_length - 1) {

			int endgap = s2_length - (c + orb);
			if (endgap > nbdf)
				return (0);

		}

	}

	//=========================================================================================

	if (((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {

		return (1);
	}

	int dist = 0;

	int flag = 0;

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{

		//cout << " S1 " << s1[c + ora] << " S2 " << s2[c + orb] << endl;

		if (endgapflag == 1) {

			//cout << " hello 2 " << endl;

			if ((c + orb >= s2_length - 1)) {

				//cout << " hello 3 : :  " << endl;

				int endgap = s1_length - (c + ora);
				if (endgap > nbdf)
					return (0);

			}

			else if (c + ora >= s1_length - 1) {

				int endgap = s2_length - (c + orb);
				if (endgap > nbdf)
					return (0);

			}

		}

		if (s1[c + ora] == s2[c + orb]) {
			if (c >= 1) {
				if (s1[c + ora] == s1[c + ora - 1]) {
					if (s1[c + ora - 1] == s2[c + orb - 1])
						counterr++;
					else
						counterr = 1;
				} else
					counterr = 1;
			} else
				counterr = 1;

			//=========================================================================================

			//=========================================================================================

			if (((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {

				return (1);

			}

			c++;

		}

		else

		{

			//cout << "Value of nbdf s1 s2 ............ dmainunmatched "
			//	<< "  Nbdf " << nbdf << " " << s1[c + ora] << " " << s2[c
			//+ orb] << " " << counterr << endl;

			if (counterr >= 2) {
				flag = 1;

			} else if (counterr == 1)
				flag = 2;
			counterr = 0;
			if ((flag == 1) && (s1[c + ora] == s1[c + ora - 1])) {
				flagi = 1;
				flagd = 0;
				flag = 0;
			} else if ((flag == 1) && (s2[c + orb] == s2[c + orb - 1])) {
				flagd = 1;
				flagi = 0;
				flag = 0;
			} else if ((flag == 2) && (s1[c + ora] == s1[c + ora - 1])) {
				flagi = 1;
				flagd = 0;
				flag = 0;
			} else if ((flag == 2) && (s2[c + orb] == s2[c + orb - 1])) {
				//if (s1[c + ora] == s1[c + ora - 1]) {

				flagd = 1;
				flagi = 0;
				flag = 0;
				/*} else {
				 flagi = 0;
				 flagd = 0;
				 flag = 0;
				 } */
			} else {
				flagi = 0;
				flagd = 0;
				flag = 0;
			}
			signed int k = (s1_length - (c + ora)) - (s2_length - (c + orb));

			//==========================================================================================

			//cout << "Valueee of  k " << flagi << " " << flagd << " " << flag
			//	<< " " << k << endl;

			if (k < 0) {

				for (int i = 1; i <= msearch; i++)

				{
					//if (c > 0) {
					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int distt = 1;

						//==========================================================================================

						/*if (flagd == 1) {
						 distt = dist;
						 }

						 else {
						 if (i > 1) {
						 distt = distt + i - 1;
						 } else
						 distt = i;
						 }
						 */

						//==========================================================================================
						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								distt = 0;
							else
								distt = i - 1;
						} else if ((i == 1) && (flagd == 1)) {

							distt = 0;
						}

						else {
							if (i > 1) {
								distt = distt + i - 1;
							} else
								distt = i;

						}

						//cout << " Value of Distance             " << distt
						//	<< endl;

						//gdist = dist;

						//==========================================================================================

						int orbb = orb + i;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (/*((abs)(kk) > nbdf - distt) ||*/(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);

							} else
								break;
						}

					}

				}

				//==========================================================================================

				if (c == 0) {

					for (int i = 1; i <= msearch; i++)

					{
						if ((c + orb + i < s2_length)

						&& (s1[c + ora + i - 1] == s2[c + orb + i]))

						{
							//cout << " i am here 1 < ++++ " << s1[c + ora + i
							//- 1] << " " << s2[c + orb + i] << endl;
							int distt = 1;

							/*if (flagd == 1) {
							 distt = dist;
							 }

							 else {
							 if (i > 1) {
							 distt = distt + i - 1;
							 } else
							 distt = i;
							 }
							 */
							//==========================================================================================
							if ((i > 1) && (flagd == 1)) {
								int flagdin = 0;
								for (int ik = 1; ik < i; ik++) {
									if (s2[c + orb] == s2[c + orb + ik])
										flagdin = 1;
									else {
										flagdin = 0;
										break;
									}
								}
								if (flagdin == 1)
									distt = 0;
								else
									distt = i - 1;
							} else if ((i == 1) && (flagd == 1)) {

								distt = 0;
							}

							else {
								if (i > 1) {
									distt = distt + i - 1;
								} else
									distt = i;

							}

							//gdist = dist;

							//==========================================================================================
							int orbb = orb + i;
							int oraa = ora + i - 1;

							signed int kk = (s1_length - c + 1 - oraa)
									- (s2_length - c + 1 - orbb);
							if (/*((abs)(kk) > nbdf - distt) ||*/(distt > nbdf))
								break;

							else {
								if (Distancemainhl(s1, s2, s1_length, s2_length,
										msearch, nbdf - distt, c + 1, oraa,
										orbb)) {

									return (1);

								} else
									break;
							}

						}

					}

				}
				//==========================================================================================

				//}

				for (int i = 1; i <= msearch; i++)

				{

					//if (c > 0) {

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{
						//cout << " i am here 2 < " << endl;
						int distt = 1;

						/*if (flagi == 1) {
						 distt = dist;
						 }

						 else {
						 if (i > 1) {
						 distt = distt + i - 1;
						 } else
						 distt = i;
						 }*/

						//===========================================================================================
						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s1[c + ora] == s1[c + ora + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								distt = 0;
							else
								distt = i - 1;

						} else if ((i == 1) && (flagi == 1)) {

							distt = 0;
						}

						else {
							if (i > 1) {
								distt = distt + i - 1;

							} else {
								distt = i;

							}
							//gdist = dist;
						}

						//===========================================================================================

						int oraa = ora + i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}

				}
				//==========================================================================================

				//==========================================================================================

				if (c == 0) {

					for (int i = 1; i <= msearch; i++) {

						if ((c + ora + i < s1_length)

						&& (s1[c + ora + i] == s2[c + orb + i - 1]))

						{
							//cout << " i am here 2 < " << endl;
							int distt = 1;

							/*if (flagi == 1) {
							 distt = dist;
							 }

							 else {
							 if (i > 1) {
							 distt = distt + i - 1;
							 } else
							 distt = i;
							 }*/

							//===========================================================================================
							if ((i > 1) && (flagi == 1)) {
								int flagdin = 0;
								for (int ik = 1; ik < i; ik++) {
									if (s1[c + ora] == s1[c + ora + ik])
										flagdin = 1;
									else {
										flagdin = 0;
										break;
									}
								}
								if (flagdin == 1)
									distt = 0;
								else
									distt = i - 1;

							} else if ((i == 1) && (flagi == 1)) {

								distt = 0;
							}

							else {
								if (i > 1) {
									distt = distt + i - 1;

								} else {
									distt = i;

								}
								//gdist = dist;
							}

							//===========================================================================================

							int oraa = ora + i;
							int orbb = orb + i - 1;

							signed int kk = (s1_length - c + 1 - oraa)
									- (s2_length - c + 1 - orbb);
							if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
								break;

							else {
								if (Distancemainhl(s1, s2, s1_length, s2_length,
										msearch, nbdf - distt, c + 1, oraa,
										orbb)) {

									return (1);
								} else
									break;
							}

						}

					}
				}

				//===========================================================================================

				//}

				for (int i = 1; i <= msearch; i++)

				{
					//cout << " i am here 3 < " << s1[c + ora + i] << " &&"
					//	<< s2[c + orb + i] << endl;
					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						//cout << " i am here 3 < " << s1[c + ora + i] << " "
						//<< s2[c + orb + i] << endl;
						int distt = 1;
						//=============================================

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								distt = 0;
							} else
								distt = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								distt = 0;
							} else
								distt = 1;

						}

						/*
						 if ((i == 1) && (flagi == 1))

						 if (s2[c + orb] == s2[c + orb + i]) {

						 distt = 0;
						 } */
						//======================================================
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = distt + i - 1;
						} else
							distt = i;
						//cout << "distt " << distt << endl;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}
				return (0);

			}

			else if (k == 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						//	cout << "Value off s1 ====== " << s1[c + ora + i]
						//		<< " " << "Value off s2 " << s2[c + orb + i]
						//	<< "  " << i << endl;

						int distt = 1;
						//=============================================
						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								distt = 0;
							} else
								distt = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								distt = 0;
							} else
								distt = 1;

						}

						/*if ((i == 1) && (flagi == 1))

						 if (s2[c + orb] == s2[c + orb + i]) {

						 distt = 0;
						 } */
						//=============================================
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = distt + i - 1;
						} else
							distt = i;

						//cout << "distt == = = " << distt << endl;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {
								// cout<<"I have returned "<<endl;
								return (1);
							} else
								break;

						}

					}

				}

				for (int i = 1; i <= msearch; i++)

				{

					//if (c > 0) {

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{

						//cout << "Value off s1 ++++" << s1[c + ora + i] << " "
						//<< "Value off s2 " << s2[c + orb] << "  "
						//<< i << endl;

						int distt = 1;

						/*if (flagi == 1) {
						 distt = dist;
						 }

						 else {
						 if (i > 1) {
						 distt = distt + i - 1;
						 } else
						 distt = i;
						 }
						 */
						//===========================================================================================
						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s1[c + ora] == s1[c + ora + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								distt = 0;
							else
								distt = i - 1;

						} else if ((i == 1) && (flagi == 1)) {

							distt = 0;
						}

						else {
							if (i > 1) {
								distt = distt + i - 1;

							} else {
								distt = i;

							}
							//gdist = dist;
						}
						//cout<<"value of distance "<<distt<<endl;
						//===========================================================================================
						int oraa = ora + i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;

						}

					}

				} //else

				//===========================================================================================
				if (c == 0) {
					for (int i = 1; i <= msearch; i++)

					{

						if ((c + ora + i < s1_length)

						&& (s1[c + ora + i] == s2[c + orb + i - 1]))

						{
							int distt = 1;

							/*if (flagi == 1) {
							 distt = dist;
							 }

							 else {
							 if (i > 1) {
							 distt = distt + i - 1;
							 } else
							 distt = i;
							 }
							 */
							//===========================================================================================
							if ((i > 1) && (flagi == 1)) {
								int flagdin = 0;
								for (int ik = 1; ik < i; ik++) {
									if (s1[c + ora] == s1[c + ora + ik])
										flagdin = 1;
									else {
										flagdin = 0;
										break;
									}
								}
								if (flagdin == 1)
									distt = 0;
								else
									distt = i - 1;

							} else if ((i == 1) && (flagi == 1)) {

								distt = 0;
							}

							else {
								if (i > 1) {
									distt = distt + i - 1;

								} else {
									distt = i;

								}
								//gdist = dist;
							}

							//===========================================================================================
							int oraa = ora + i;
							int orbb = orb + i - 1;

							signed int kk = (s1_length - c + 1 - oraa)
									- (s2_length - c + 1 - orbb);
							if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
								break;

							else {
								if (Distancemainhl(s1, s2, s1_length, s2_length,
										msearch, nbdf - distt, c + 1, oraa,
										orbb)) {

									return (1);
								} else
									break;

							}

						}

					}
					//==========================================================================================

				}

				for (int i = 1; i <= msearch; i++) {

					//if (c > 0) {
					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{
						//cout << "Value off s1 ****" << s1[c + ora] << " "
						//<< "Value off s2 " << s2[c + orb + i] << "  "
						//<< i << endl;

						int distt = 1;

						/*if (flagd == 1) {
						 distt = dist;
						 }

						 else {
						 if (i > 1) {
						 distt = distt + i - 1;
						 } else
						 distt = i;
						 }*/
						//==========================================================================================
						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								distt = 0;
							else
								distt = i - 1;
						} else if ((i == 1) && (flagd == 1)) {

							distt = 0;
						}

						else {
							if (i > 1) {
								distt = distt + i - 1;
							} else
								distt = i;

						}
						//==========================================================================================
						int orbb = orb + i;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {
								//cout<<"Hello  world "<<endl;
								return (1);
							} else
								break;

						}

					}
				} //else
				  //==========================================================================================
				if (c == 0) {
					for (int i = 1; i <= msearch; i++)

					{

						if ((c + orb + i < s2_length)

						&& (s1[c + ora + i - 1] == s2[c + orb + i]))

						{
							int distt = 1;

							/*if (flagd == 1) {
							 distt = dist;
							 }

							 else {
							 if (i > 1) {
							 distt = distt + i - 1;
							 } else
							 distt = i;
							 }*/
							//==========================================================================================
							if ((i > 1) && (flagd == 1)) {
								int flagdin = 0;
								for (int ik = 1; ik < i; ik++) {
									if (s2[c + orb] == s2[c + orb + ik])
										flagdin = 1;
									else {
										flagdin = 0;
										break;
									}
								}
								if (flagdin == 1)
									distt = 0;
								else
									distt = i - 1;
							} else if ((i == 1) && (flagd == 1)) {

								distt = 0;
							}

							else {
								if (i > 1) {
									distt = distt + i - 1;
								} else
									distt = i;

							}
							//==========================================================================================
							int orbb = orb + i;
							int oraa = ora + i - 1;

							signed int kk = (s1_length - c + 1 - oraa)
									- (s2_length - c + 1 - orbb);
							if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
								break;

							else {
								if (Distancemainhl(s1, s2, s1_length, s2_length,
										msearch, nbdf - distt, c + 1, oraa,
										orbb)) {

									return (1);
								} else
									break;

							}

						}

					}
					//==========================================================================================
				}
				return (0);

			}

			else if (k > 0) {

				for (int i = 1; i <= msearch; i++)

				{

					//if (c > 0) {
					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{
						//cout << " i am here 1 " << endl;
						//cout << " i am here 1  > > > " << s1[c + ora + i]
						//	<< " " << s2[c + orb] << " c " << c << " ora "
						//<< ora << " orb " << orb << endl;

						int distt = 1;

						/*if (flagi == 1) {
						 distt = dist;
						 }

						 else {
						 if (i > 1) {
						 distt = distt + i - 1;
						 } else
						 distt = i;
						 }*/

						//===========================================================================================
						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s1[c + ora] == s1[c + ora + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								distt = 0;
							else
								distt = i - 1;

							//cout << "value of distance lll^ " << i << " "
							//	<< distt << endl;

						} else if ((i == 1) && (flagi == 1)) {

							distt = 0;

						}

						else {
							if (i > 1) {
								distt = distt + i - 1;
								//cout << "hhhhh " << endl;
							} else {
								distt = i;

							}
							//gdist = dist;
						}

						//cout << "value of distance llll " << i << " " << distt
						//	<< endl;
						//cout << "value of nbdf  " << nbdf << endl;

						//===========================================================================================

						int oraa = ora + i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}
				}
				//==========================================================================================
				//else
				if (c == 0) {
					for (int i = 1; i <= msearch; i++)

					{

						if ((c + ora + i < s1_length)

						&& (s1[c + ora + i] == s2[c + orb + i - 1]))

						{
							//cout << " i am here 1 " << endl;

							int distt = 1;

							/*if (flagi == 1) {
							 distt = dist;
							 }

							 else {
							 if (i > 1) {
							 distt = distt + i - 1;
							 } else
							 distt = i;
							 }*/

							//===========================================================================================
							if ((i > 1) && (flagi == 1)) {
								int flagdin = 0;
								for (int ik = 1; ik < i; ik++) {
									if (s1[c + ora] == s1[c + ora + ik])
										flagdin = 1;
									else {
										flagdin = 0;
										break;
									}
								}
								if (flagdin == 1)
									distt = 0;
								else
									distt = i - 1;

							} else if ((i == 1) && (flagi == 1)) {

								distt = 0;
							}

							else {
								if (i > 1) {
									distt = distt + i - 1;

								} else {
									distt = i;

								}
								//gdist = dist;
							}

							//===========================================================================================

							int oraa = ora + i;
							int orbb = orb + i - 1;

							signed int kk = (s1_length - c + 1 - oraa)
									- (s2_length - c + 1 - orbb);
							if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
								break;

							else {
								if (Distancemainhl(s1, s2, s1_length, s2_length,
										msearch, nbdf - distt, c + 1, oraa,
										orbb)) {

									return (1);
								} else
									break;
							}

						}
					}

					//==========================================================================================
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						/*cout << " i am here 2 " << s1[c + ora + i] << " "
						 << s2[c + orb + i] << " " << i << " " << nbdf
						 << endl;*/

						int distt = 1;
						//=====================================================
						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								distt = 0;
							} else
								distt = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								distt = 0;
							} else
								distt = 1;

						}

						/*if ((i == 1) && (flagi == 1))

						 if (s2[c + orb] == s2[c + orb + i]) {

						 distt = 0;
						 } */
						//=====================================================
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = distt + i - 1;
						} else
							distt = i;

						//cout << "distt " << distt << endl;

						//}
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (/*((abs)(kk) > nbdf - distt) ||*/(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}

				for (int i = 1; i <= msearch; i++) {

					//if (c > 0) {

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{
						//cout << " i am here 3 " << endl;
						int distt = 1;
						//==========================================================================================

						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								distt = 0;
							else
								distt = i - 1;
						} else if ((i == 1) && (flagd == 1)) {

							distt = 0;
						}

						else {
							if (i > 1) {
								distt = distt + i - 1;
							} else
								distt = i;

						}

						//==========================================================================================

						/*if (flagd == 1) {
						 distt = dist;
						 }

						 else {
						 if (i > 1) {
						 distt = distt + i - 1;
						 } else
						 distt = i;
						 }*/

						int orbb = orb + i;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
							break;

						else {
							if (Distancemainhl(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);
							} else
								break;
						}

					}
				}
				//=======================================================================

				//else
				if (c == 0) {
					for (int i = 1; i <= msearch; i++)

					{

						if ((c + orb + i < s2_length)

						&& (s1[c + ora + i - 1] == s2[c + orb + i]))

						{
							//cout << " i am here 3 " << endl;
							int distt = 1;
							//==========================================================================================

							if ((i > 1) && (flagd == 1)) {
								int flagdin = 0;
								for (int ik = 1; ik < i; ik++) {
									if (s2[c + orb] == s2[c + orb + ik])
										flagdin = 1;
									else {
										flagdin = 0;
										break;
									}
								}
								if (flagdin == 1)
									distt = 0;
								else
									distt = i - 1;
							} else if ((i == 1) && (flagd == 1)) {

								distt = 0;
							}

							else {
								if (i > 1) {
									distt = distt + i - 1;
								} else
									distt = i;

							}

							//==========================================================================================

							/*if (flagd == 1) {
							 distt = dist;
							 }

							 else {
							 if (i > 1) {
							 distt = distt + i - 1;
							 } else
							 distt = i;
							 }*/

							int orbb = orb + i;
							int oraa = ora + i - 1;

							signed int kk = (s1_length - c + 1 - oraa)
									- (s2_length - c + 1 - orbb);
							if (/*((abs)(kk) > nbdf - distt) || */(distt > nbdf))
								break;

							else {
								if (Distancemainhl(s1, s2, s1_length, s2_length,
										msearch, nbdf - distt, c + 1, oraa,
										orbb)) {

									return (1);
								} else
									break;
							}

						}

					}

					//=======================================================================
				}
				return (0);

			}

		}

	}

}
bool Distanceonehl(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb)

		{

	/*cout << "I am here ++++++++++++++ => nbdf " << nbdf << " " << c << " "
	 << ora << " " << orb << " " << endl; */
	if (endgapflag == 0) {
		if ((c + ora == s1_length) || (c + orb == s2_length))
			return (1);
	} else if (endgapflag == 1) {
		if ((c + ora == s1_length) && (c + orb == s2_length))
			return (1);
	}

	int s11_length = s1_length - c - ora;
	int s22_length = s2_length - c - orb;

	int k = (s11_length - s22_length);

	//cout << "Value of k " << k << endl;

	//=========================================================================================
	if (c + orb >= s2_length - 1) {

		if (endgapflag == 1) {
			int endgap = s1_length - (c + ora);
			if (endgap > nbdf)
				return (0);
		} else
			return (1);

	} else if (c + ora >= s1_length - 1) {

		if (endgapflag == 1) {
			int endgap = s2_length - (c + orb);
			if (endgap > nbdf)
				return (0);
		} else
			return (1);
	}

	//=========================================================================================

	int common = 0;

	int counterr;

	int flagi = 0;
	int flagd = 0;

	int flag;

	counterr = 1;

	int gdist = 0;

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{

		//cout<<"Hello   "<<endl;

		if (c + orb >= s2_length - 1) {
			if (endgapflag == 1) {
				int endgap = s1_length - (c + ora);
				if (endgap > nbdf)
					return (0);

			} else
				return (1);

		}

		else if (c + ora >= s1_length - 1) {
			if (endgapflag == 1) {

				int endgap = s2_length - (c + orb);
				if (endgap > nbdf)
					return (0);

			} else
				return (1);

		}

		//else
		//===============================================================================================

		//===============================================================================================

		//}

		/*cout << " Start <<>> S1 " << s1[c + ora] << " S2 " << s2[c + orb]
		 << "  ** " << nbdf << endl; */

		if (s1[c + ora] == s2[c + orb]) {
			if (s1[c + ora] == s1[c + ora - 1]) {
				if (s1[c + ora - 1] == s2[c + orb - 1])
					counterr++;
				else
					counterr = 1;
			} else
				counterr = 1;

			//=======================================================================

			if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
				return (1);

			//cout << " S1 ==== " << s1[c + ora] << " S2 ==== " << s2[c + orb] << "  "
			//<< counterr << endl;

			//==========================================================================================

			if ((c + ora < s1_length) && (c + orb == s2_length - 1)) {

				//=============================================================================
				if (endgapflag == 1) {
					int endgap = s1_length - (c + ora);
					if (endgap > nbdf)
						return (0);
				}

				//=============================================================================

				int flagint = 0;

				if (counterr >= 2)
					flagint = 1;
				else
					flagint = 0;

				if (flagint == 0)
					return (0);
				else if (flagint == 1) {
					for (int ik = c + ora + 1; ik < s1_length; ik++)
						if (s1[ik] == s2[c + orb])
							flagint = 1;
						else {
							gdist++;
							flagint = 0;
							break;
						}
				}

				if ((flagint == 1) && (gdist <= 1))
					return (1);
				else
					return (0);

			} else if ((c + ora == s1_length - 1) && (c + orb < s2_length)) {
				//=============================================================================
				if (endgapflag == 1) {

					int endgap = s2_length - (c + orb);
					if (endgap > nbdf)
						return (0);

				}
				//=============================================================================

				int flagint = 0;

				if (counterr >= 2)
					flagint = 1;
				else
					flagint = 0;

				if (flagint == 0)
					return (0);

				if (flagint == 1) {
					for (int ik = c + orb + 1; ik < s2_length; ik++)
						if (s1[c + ora] == s2[ik])
							flagint = 1;
						else {
							gdist++;
							flagint = 0;
							break;
						}
				}
				if ((flagint == 1) && (gdist <= 1))
					return (1);
				else
					return (0);
			}
			//==========================================================================================

			common++;
			c++;

		} else

		{

			/*cout << " S1 >>>> ==== <<<< " << s1[c + ora] << " S2 ==== > "
			 << s2[c + orb] << "  " << counterr << " " << " " << k
			 << " " << nbdf << endl; */

			// cout << " S1 >>>> ==== <<<< " << s1[c + ora] << " S2 > ==== " << s2[c + orb +1] << "  "
			//<< counterr << endl;
			if (counterr >= 2) {
				flag = 1;

			} else if (counterr == 1)
				flag = 2;
			counterr = 0;
			if ((flag == 1) && (s1[c + ora] == s1[c + ora - 1])) {
				flagi = 1;
				flagd = 0;
				flag = 0;
			} else if ((flag == 1) && (s2[c + orb] == s2[c + orb - 1])) {
				flagd = 1;
				flagi = 0;
				flag = 0;
			} else if ((flag == 2) && (s1[c + ora] == s1[c + ora - 1])) {
				flagi = 1;
				flagd = 0;
				flag = 0;
			} else if ((flag == 2) && (s2[c + orb] == s2[c + orb - 1])) {
				//if (s1[c + ora] == s1[c + ora - 1]) {

				flagd = 1;
				flagi = 0;
				flag = 0;
				/*} else {
				 flagi = 0;
				 flagd = 0;
				 flag = 0;
				 }*/
			} else {
				flagi = 0;
				flagd = 0;
				flag = 0;
			}

			//cout << " Flags and k FlaPOOOO " << flag << " S2 ==== > " << flagi
			//	<< "  " << flagd << endl;

			if (k == 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						//cout << "Value of i gdist " << i << " " << gdist
						//	<< endl;

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist += 0;
							} else
								gdist += 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist += 0;
							} else
								gdist += 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist += dist;

						}

						/*float dist = 1;

						 if (i > 1) {
						 dist = dist + i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;*/

						//cout << "Value of Gdisuiiiot " << i << " " << gdist
						//<< " " << s1[c + ora + i] << " " << s2[c + orb
						//+ i] << endl;
						if ((gdist > nbdf)
								|| (s1[c + ora + i + 1] != s2[c + orb + i + 1])) {

							goto level11;
						}

						else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							orb = orb + i;
							common++;
							goto level1;

						}

					}
					if (i == msearch)
						//gdist++;

						goto level11;

				}

				//=========================================================================================
				level11: for (int i = 1; i <= msearch; i++)

				{

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						//cout << "Value of Gdistii "<<" "<<i<<" " << gdist << " "
						//<< s1[c + ora ] << " " << s2[c + orb + i]
						//<< endl;

						float dist = 0;

						/*if (flagd == 1) {
						 dist = 0;
						 }

						 else {
						 if (i > 1) {
						 dist = i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;
						 }
						 if (gdist > nbdf) {
						 return (0);
						 }*/
						//=========================================================================================
						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist++;
						} else if ((i == 1) && (flagd == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;

						}
						gdist = dist;
						if (gdist > nbdf) {
							goto level12;
						}

						//=========================================================================================
						else {
							nbdf = nbdf - gdist;
							orb = orb + i;

							common++;

							goto level1;
						}

					}

					if (i == msearch)
						//gdist++;

						goto level12;

				}

				//=========================================================================================
				level12: for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{

						//cout << "Value of Gdist ikk" <<" "<<i<<" " << gdist << " "
						//<< s1[c + ora + i] << " " << s2[c + orb]
						//<< endl;

						float dist = 0;

						//==========================================================================================
						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s1[c + ora] == s1[c + ora + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist++;

						} else if ((i == 1) && (flagi == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;

							} else {
								dist = i + gdist;

							}

						}
						gdist = dist;
						//cout << "gdist " << gdist << endl;

						if (gdist > nbdf) {
							return (0);
						}

						//==========================================================================================

						/*if (flagi == 1) {
						 dist = 0;
						 }

						 else {
						 if (i > 1) {
						 dist = i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;
						 }
						 if (gdist > nbdf) {
						 return (0);
						 }
						 */
						else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							common++;

							goto level1;
						}

					}
					if (i == msearch)
						gdist++;
					if (gdist > nbdf) {
						return (0);
					} else {
						nbdf = nbdf - gdist;
						goto level1;
					}

				}
				level1: c++;

				//=========================================================================================
			}

			//===========================================================================================

			else if (k < 0) {

				//===============================================================

				for (int i = 1; i <= msearch; i++)

				{
					//cout<<"Hello "<<endl;

					//cout << "Value of gdist nbdf " << nbdf << " " << ora << " "
					//	<< orb << " " << gdist << endl;

					//cout << "Value of Gdist Poinnq == i<< " << i << " "
					//	<< gdist << " " << s1[c + ora + i] << " " << s2[c
					//+ orb] << endl;

					//cout << " === s dfsssd " << s1[c + ora] << " " << s2[c
					//	+ orb] << endl;

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist = 0;
							} else
								gdist = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist += 0;
							} else
								gdist += 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist += dist;
						}

						/*float dist = 1;

						 if (i > 1) {
						 dist = dist + i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist; */

						//cout << "Value of Gdist wÃ¹Ã¹Ã¹ " << gdist << " " << s1[c
						//	+ ora + i] << " " << s2[c + orb + i] << gdist
						//<< nbdf << endl;
						if ((gdist > nbdf)
								|| (s1[c + ora + i + 1] != s2[c + orb + i + 1])) {
							goto level50;
							//return (0);
						} else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							orb = orb + i;
							common++;
							goto level2;
						}

					}

					//cout<<i<<" "<<msearch<<endl;

					if ((gdist > nbdf) || (i == 1)) {
						//break;
						goto level50;
						//return (0);
					} else {
						nbdf = nbdf - gdist;
						goto level2;
					}

					if (i == msearch) {
						//gdist++;
						goto level50;
					}

				}

				//================================================================

				level50: for (int i = 1; i <= msearch; i++)

				{
					//cout << "====== Wowowwo " << endl;

					//cout << " === s >>> " << s1[c + ora] << " " << s2[c + orb
					//	+ i] << "=== " << s2[c + orb] << " " << s2[c + orb
					//- 1] << flagd << endl;
					//=============================================================
					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						float dist = 0;

						//cout << "Value of Flagd " << flagd << endl;

						/*if (flagd == 1) {
						 dist = 0;
						 }

						 else {
						 if (i > 1) {
						 dist = i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;
						 }
						 if (gdist > nbdf) {
						 return (0);
						 }*/
						//=========================================================================================
						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist = i - 1;
						} else if ((i == 1) && (flagd == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;

						}
						gdist = dist;

						//cout << "GGGDist " << i << " " << dist << endl;

						if (gdist > nbdf) {
							goto level21;
						}

						//=========================================================================================
						else {
							nbdf = nbdf - gdist;
							orb = orb + i;

							common++;

							goto level2;
						}

					}

					if (i == msearch) {
						//gdist++;
						goto level21;
					}

				}
				//=====================================================================================

				level21: for (int i = 1; i <= msearch; i++)

				{

					//cout << "====== Wowowwo 2 " << endl;

					//cout << " === s >>> " << s1[c + ora + i] << " " << s2[c
					//	+ orb] << endl;

					if ((c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{
						//cout << gdist << endl;

						float dist = 0;

						//==========================================================================================
						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s1[c + ora] == s1[c + ora + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist = i - 1;

						} else if ((i == 1) && (flagi == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;

							} else {
								dist = i + gdist;

							}

						}
						//==========================================
						gdist = dist;
						//===================================
						//	cout << "gdist " << gdist << endl;

						if (gdist > nbdf) {
							goto level22;
						}

						//==========================================================================================

						/*if (flagi == 1) {
						 dist = 0;
						 }

						 else {
						 if (i > 1) {
						 dist = i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;
						 }
						 if (gdist > nbdf) {
						 return (0);
						 }
						 */
						else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							common++;

							goto level2;
						}

					}
					if (i == msearch) {
						//gdist++;
						goto level22;
					}

				}
				//==========================================================================================

				level22: for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist = 0;
							} else
								gdist = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist += 0;
							} else
								gdist += 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist += dist;
						}

						/*float dist = 1;

						 if (i > 1) {
						 dist = dist + i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist; */

						//cout << "Value of Gdist " << gdist << " "
						//<< s1[c + ora + i] << " " << s2[c + orb + i]
						//<< endl;
						if (gdist > nbdf) {

							return (0);
						} else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							orb = orb + i;
							common++;
							goto level2;
						}

					}
					if (i == msearch)
						gdist++;
					if (gdist > nbdf) {

						return (0);
					} else {
						nbdf = nbdf - gdist;
						goto level2;
					}

				}

				level2: c++;
				//===========================================================================================
			}

			//==========================================================================================
			else if (k > 0) {

				for (int i = 1; i <= msearch; i++)

				{

					//cout << "Value of gdist nbdf " << nbdf << " " << ora << " "
					//	<< orb << " " << gdist << endl;

					//cout << "Value of Gdist Poinnq == i " << i << " " << gdist
					//	<< " " << s1[c + ora + i + 1] << " " << s2[c + orb]
					//<< endl;

					//cout << " === s " << s1[c + ora] << " " << s2[c + orb]
					//	<< endl;
					/*cout << " === s " << s1[c + ora - 1] << " " << s2[c + orb
					 - 1] << endl;
					 cout << " === s " << s1[c + ora - 2] << " " << s2[c + orb
					 - 2] << endl;*/

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						//==================================================

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist = 0;
							} else
								gdist = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist += 0;
							} else
								gdist += 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist += dist;
						}

						//cout << "Value of Gdist Ã Ã Ã  " << gdist << " " << s1[c
						//	+ ora + i] << " " << s2[c + orb + i] << " == "
						//<< nbdf << endl;

						if ((gdist > nbdf)
								|| (s1[c + ora + i + 1] != s2[c + orb + i + 1])) {
							goto level30;
							//return (0);
						} else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							orb = orb + i;
							common++;
							goto level3;
						}

					}
					if (i == msearch)
						gdist++;
					if ((gdist > nbdf) || (i == 1)) {
						goto level30;
						//return (0);
					} else
						goto level3;

					if (i == msearch) {
						gdist++;
						//cout << "======================== " << endl;
						goto level30;
					}

				}

				level30: for (int i = 1; i <= msearch; i++)

				{

					//cout << "Hello  ____ >><<< " << s1[c + ora + 2] << s2[c
					//	+ orb] << endl;

					if ((c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{
						float dist = 0;

						//cout << "Hello  ____ " << endl;

						//==========================================================================================
						if ((i > 1) && (flagi == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s1[c + ora] == s1[c + ora + ik])
								//==================================
										{

									flagdin = 1;
								}
								//=====================================

								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist = i - 1;

						} else if ((i == 1) && (flagi == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;

							} else {
								dist = i + gdist;

							}

						}
						gdist = dist;

						//cout << "gdist nbdf ::: " << gdist << " " << nbdf
						//	<< endl;

						if (gdist > nbdf) {
							goto level31;
						}

						//==========================================================================================

						/*if (flagi == 1) {
						 dist = 0;
						 }

						 else {
						 if (i > 1) {
						 dist = i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;
						 }
						 if (gdist > nbdf) {
						 return (0);
						 }
						 */
						else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							common++;
							//cout << "gdist nbdf ::: >>>" << gdist << " "
							//	<< nbdf << endl;

							//cout << s1[c + ora] << " " << s2[c + orb] << endl;

							goto level3;
						}

					}
					if (i == msearch) {
						gdist++;
						//cout << "======================== " << endl;
						goto level31;
					}

				}

				//=====================================================================================
				level31: for (int i = 1; i <= msearch; i++)

				{
					//cout << "========================l:l " << endl;
					//cout << "Hello  ____ >>>>> " << s1[c + ora] << " " << s2[c
					//	+ orb + i] << " " << nbdf << endl;

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{
						float dist = 0;

						//cout << "========================l:l " << ora << orb
						//	<< endl;

						/*if (flagd == 1) {
						 dist = 0;
						 }

						 else {
						 if (i > 1) {
						 dist = i - 1 + gdist;
						 } else
						 dist = i + gdist;
						 gdist = dist;
						 }
						 if (gdist > nbdf) {
						 return (0);
						 }*/
						//=========================================================================================
						if ((i > 1) && (flagd == 1)) {
							int flagdin = 0;
							for (int ik = 1; ik < i; ik++) {
								if (s2[c + orb] == s2[c + orb + ik])
									flagdin = 1;
								else {
									flagdin = 0;
									break;
								}
							}
							if (flagdin == 1)
								dist = 0;
							else
								dist = i - 1;
						} else if ((i == 1) && (flagd == 1)) {

							dist = 0;
						}

						else {
							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i;

						}
						gdist = dist;

						//cout << "Value of gdist i " << i << ora << orb << gdist
						//	<< endl;
						if (gdist > nbdf) {

							/*cout << "Value of gdist " << ora << orb << gdist
							 << endl;*/

							goto level32;
						}

						//=========================================================================================
						else {
							nbdf = nbdf - gdist;
							orb = orb + i;

							common++;

							goto level3;
						}

					}

					if (i == msearch) {
						gdist++;
						//cout << "========================  Ã© " << endl;
						goto level32;
					}

				}

				//=========================================================================================
				level32: for (int i = 1; i <= msearch; i++)

				{

					/*cout << "Value of gdist " << ora << orb << gdist << endl;

					 cout << "Value of Gdist == i " << i << " " << gdist << " "
					 << s1[c + ora + i] << " " << s2[c + orb + i]
					 << endl;

					 cout << " === s " << s1[c + ora] << " " << s2[c + orb]
					 << endl;
					 cout << " === s " << s1[c + ora - 1] << " " << s2[c + orb
					 - 1] << endl;
					 cout << " === s " << s1[c + ora - 2] << " " << s2[c + orb
					 - 2] << endl;*/

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{

						//==================================================

						/*if ((flagi == 1)&&(s2[c + orb] == s2[c + orb + i])) {
						 //if (gdist >= 2)
						 gdist = 0;
						 } else*/

						if ((i == 1) && (flagi == 1)) {

							if (s2[c + orb] == s2[c + orb + i]) {

								gdist = 0;
							} else
								gdist = 1;
						} else if ((i == 1) && (flagd == 1)) {

							if (s1[c + ora] == s1[c + ora + i]) {

								gdist = 0;
							} else
								gdist = 1;

						} else {
							float dist = 1;

							//==================================================

							if (i > 1) {
								dist = dist + i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist = dist;
						}

						//cout << "Value of Gdist Ã Ã Ã  " << gdist << " " << s1[c
						//	+ ora + i] << " " << s2[c + orb + i] << endl;

						if (gdist > nbdf) {

							return (0);
						} else {
							nbdf = nbdf - gdist;
							ora = ora + i;
							orb = orb + i;
							common++;
							goto level3;
						}

					}
					if (i == msearch)
						gdist++;
					if (gdist > nbdf) {

						return (0);
					} else
						goto level3;

				}

				//=========================================================================================
				level3: c++;
			}
			//c++;

		}

	}

	return (1);

}

//==================================================================================================================

//========================================== PointMutation+Homologues+HH=>example( CG, CC, AA, AT )======================

bool Distancehh(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int seed_id_size, int test_id_size, int msearch,
		int nbdf)

		{

	int k = (s1_length - seed_id_size) - (s2_length - test_id_size);

	int c = 0;

	int ora = seed_id_size;

	int orb = test_id_size;

	int common = 0;
	int flag;
	int counterr = 0;
	float gdist = 0;

	if (nbdf == 1) {
		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{
			flag = 0;

			if (s1[c + ora] == s2[c + orb]) {
				//=========================================================================================
				if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
					return (1);

				else if ((c + ora < s1_length) && (c + orb == s2_length - 1)) {

					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);
					else if (flagint == 1) {
						for (int ik = c + ora + 1; ik < s1_length; ik++)
							if (s1[c + ora + ik] == s2[c + orb])
								flagint = 1;
							else {
								gdist++;
								flagint = 0;
								break;
							}
					}

					if ((flagint == 1) && (gdist <= 1))
						return (1);
					else
						return (0);

				} else if ((c + ora == s1_length - 1)
						&& (c + orb < s2_length)) {
					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);

					if (flagint == 1) {
						for (int ik = c + orb + 1; ik < s2_length; ik++)
							if (s1[c + ora] == s2[c + orb + ik])
								flagint = 1;
							else {
								gdist++;
								flagint = 0;
								break;
							}
					}
					if ((flagint == 1) && (gdist <= 1))
						return (1);
					else
						return (0);
				}

				//==========================================================================================

				counterr++;
				common++;
				c++;

			} else

			{

				if (counterr >= 2) {
					flag = 1;
				}
				counterr = 0;

				for (int i = 1; i <= msearch; i++)

				{

					if ((k == 0) && (c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						counterr++;
						float dist = 0;
						ora = ora + i;
						orb = orb + i;
						if (i > 1) {
							dist = i - 1 + gdist;

						} else
							dist = i + gdist;
						gdist = dist;

						if (gdist > nbdf) {
							return (0);
						}

						common++;

						break;

					}

					else if ((k < 0) && (c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{
						counterr++;
						float dist = 0;

						if (flag == 1) {
							dist = 0;
						}

						else {
							if (i > 1) {
								dist = i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist = dist;
						}
						if (gdist > nbdf) {
							return (0);
						}

						//====================================================================================

						orb = orb + i;

						common++;

						break;

					}

					else if ((k > 0) && (c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{
						counterr++;
						float dist = 0;

						if (flag == 1) {
							dist = 0;
						}

						else {
							if (i > 1) {
								dist = i - 1 + gdist;
							} else
								dist = i + gdist;
							gdist = dist;
						}

						if (gdist > nbdf) {
							return (0);
						}

						ora = ora + i;

						common++;

						break;

					}
					if (i == msearch)
						gdist++;

				}
				c++;

			}

		}

		return (0);

	}
	//====================================[ NBDF==0 ]================================================

	else if (nbdf == 0) {
		while ((c + ora < s1_length)

		&& (c + orb < s2_length))

		{
			flag = 0;

			if (s1[c + ora] == s2[c + orb]) {

				//==========================================================================================

				if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
					return (1);

				else if ((c + ora < s1_length) && (c + orb == s2_length - 1)) {

					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);
					if (flagint == 1) {
						for (int ik = c + ora + 1; ik < s1_length; ik++)
							if (s1[c + ora + ik] == s2[c + orb])
								flagint = 1;
							else
								flagint = 0;
					}

					if (flagint == 1)
						return (1);
					else
						return (0);

				} else if ((c + ora == s1_length - 1)
						&& (c + orb < s2_length)) {
					int flagint = 0;

					if (counterr >= 2)
						flagint = 1;
					else
						flagint = 0;

					if (flagint == 0)
						return (0);

					if (flagint == 1) {
						for (int ik = c + orb + 1; ik < s2_length; ik++)
							if (s1[c + ora] == s2[c + orb + ik])
								flagint = 1;
							else
								flagint = 0;
					}
					if (flagint == 1)
						return (1);
					else if (flagint == 0)
						return (0);
				}

				//==========================================================================================

				counterr++;
				common++;
				c++;

			} else

			{

				if (counterr >= 2) {
					flag = 1;
				}
				counterr = 0;

				for (int i = 1; i <= msearch; i++)

				{

					if ((k == 0) && (c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						counterr++;
						float dist = 0;
						ora = ora + i;
						orb = orb + i;
						if (i > 1) {
							dist = i - 1;
						} else
							dist = i;

						if (dist > nbdf) {
							return (0);
						}

						common++;

						break;

					}

					else if ((k < 0) && (c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{
						counterr++;
						float dist = 0;

						if (flag == 1) {
							dist = 0;
						}

						else {
							if (i > 1) {
								dist = i - 1;
							} else
								dist = i;
						}
						if (dist > nbdf) {
							return (0);
						}

						//====================================================================================

						orb = orb + i;

						common++;

						break;

					}

					else if ((k > 0) && (c + ora + i < s1_length)
							&& (s1[c + ora + i] == s2[c + orb]))

							{
						counterr++;
						float dist = 0;

						if (flag == 1) {
							dist = 0;
						}

						else {
							if (i > 1) {
								dist = i - 1;
							} else
								dist = i;
						}

						if (dist > nbdf) {
							return (0);
						}

						ora = ora + i;

						common++;

						break;

					}

					if (i == msearch)
						return (0);
				}
				c++;

			}

		}

		return (0);

	}

	//=========================================================================================

	else {
		globalflag = 1;
		if (Distancemainhh(s1, s2, s1_length, s2_length, msearch, nbdf, 0,
				seed_id_size, test_id_size))
			return (1);
		else
			return (0);
	}

}

bool Distancemainhh(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb) {

	int counterr;
	if (globalflag == 1) {
		counterr = 0;
		globalflag = 0;
	} else
		counterr = 1;

	if (nbdf <= 1) {

		if (Distanceonehh(s1, s2, s1_length, s2_length, msearch, nbdf, c, ora,
				orb))
			return (1);
		else
			return (0);
	}

	signed int kk = (s1_length - c - ora) - (s2_length - c - orb);

	if (((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {

		return (1);
	}

	int dist = 0;

	int flag;

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{
		flag = 0;

		if (s1[c + ora] == s2[c + orb]) {
			counterr++;

			if (((c + ora == s1_length - 1) || (c + orb == s2_length - 1))) {

				return (1);

			}

			c++;

		}

		else

		{
			if (counterr >= 2) {
				flag = 1;
			}
			counterr = 0;

			signed int k = (s1_length - (c + ora)) - (s2_length - (c + orb));

			if (k < 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int distt = 0;

						if (flag == 1) {
							distt = dist;
						}

						else {
							if (i > 1) {
								distt = i - 1;
							} else
								distt = i;
						}

						int orbb = orb + i;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);

							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{
						int distt = 0;

						if (flag == 1) {
							distt = dist;
						}

						else {
							if (i > 1) {
								distt = i - 1;
							} else
								distt = i;
						}

						int oraa = ora + i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						int distt = 0;
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = i - 1;
						} else
							distt = i;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}
				return (0);

			}

			else if (k == 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						int distt = 0;
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = i - 1;
						} else
							distt = i;
						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;

						}

					}

				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{
						int distt = 0;

						if (flag == 1) {
							distt = dist;
						}

						else {
							if (i > 1) {
								distt = i - 1;
							} else
								distt = i;
						}

						int oraa = ora + i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;

						}

					}
				}

				for (int i = 1; i <= msearch; i++) {

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{
						int distt = 0;

						if (flag == 1) {
							distt = dist;
						}

						else {
							if (i > 1) {
								distt = i - 1;
							} else
								distt = i;
						}

						int orbb = orb + i;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);
							} else
								break;

						}

					}
				}
				return (0);

			}

			else if (k > 0) {

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb]))

					{
						int distt = 0;

						if (flag == 1) {
							distt = dist;
						}

						else {
							if (i > 1) {
								distt = i - 1;
							} else
								distt = i;
						}

						int oraa = ora + i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orb)) {

								return (1);
							} else
								break;
						}

					}
				}

				for (int i = 1; i <= msearch; i++)

				{

					if ((c + ora + i < s1_length)

					&& (s1[c + ora + i] == s2[c + orb + i]))

					{
						int distt = 0;
						int oraa = ora + i;
						int orbb = orb + i;
						if (i > 1) {
							distt = i - 1;
						} else
							distt = i;

						signed int kk = (s1_length - c + 1 - oraa)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, oraa, orbb)) {

								return (1);
							} else
								break;
						}

					}

				}

				for (int i = 1; i <= msearch; i++) {

					if ((c + orb + i < s2_length)

					&& (s1[c + ora] == s2[c + orb + i]))

					{

						int distt = 0;

						if (flag == 1) {
							distt = dist;
						}

						else {
							if (i > 1) {
								distt = i - 1;
							} else
								distt = i;
						}

						int orbb = orb + i;

						signed int kk = (s1_length - c + 1 - ora)
								- (s2_length - c + 1 - orbb);
						if (((abs)(kk) > nbdf - distt) || (distt > nbdf))
							break;

						else {
							if (Distancemainhh(s1, s2, s1_length, s2_length,
									msearch, nbdf - distt, c + 1, ora, orbb)) {

								return (1);
							} else
								break;
						}

					}
				}
				return (0);

			}

		}

	}

}

bool Distanceonehh(vector<char>& s1, vector<char>& s2, int s1_length,
		int s2_length, int msearch, int nbdf, int c, int ora, int orb)

		{

	int s11_length = s1_length - c - ora;
	int s22_length = s2_length - c - orb;

	int k = (s11_length - s22_length);

	int common = 0;

	int counterr;

	int flag;

	counterr = 1;

	float gdist = 0;

	while ((c + ora < s1_length)

	&& (c + orb < s2_length))

	{
		flag = 0;

		if (s1[c + ora] == s2[c + orb]) {
			counterr++;

			if ((c == s1_length - ora - 1) || (c == s2_length - orb - 1)) {
				return (1);

			}
			//==========================================================================================
			if ((c + ora == s1_length - 1) && (c + orb == s2_length - 1))
				return (1);

			else if ((c + ora < s1_length) && (c + orb == s2_length - 1)) {

				int flagint = 0;

				if (counterr >= 2)
					flagint = 1;
				else
					flagint = 0;

				if (flagint == 0)
					return (0);
				else if (flagint == 1) {
					for (int ik = c + ora + 1; ik < s1_length; ik++)
						if (s1[c + ora + ik] == s2[c + orb])
							flagint = 1;
						else {
							gdist++;
							flagint = 0;
							break;
						}
				}

				if ((flagint == 1) && (gdist <= 1))
					return (1);
				else
					return (0);

			} else if ((c + ora == s1_length - 1) && (c + orb < s2_length)) {
				int flagint = 0;

				if (counterr >= 2)
					flagint = 1;
				else
					flagint = 0;

				if (flagint == 0)
					return (0);

				if (flagint == 1) {
					for (int ik = c + orb + 1; ik < s2_length; ik++)
						if (s1[c + ora] == s2[c + orb + ik])
							flagint = 1;
						else {
							gdist++;
							flagint = 0;
							break;
						}
				}
				if ((flagint == 1) && (gdist <= 1))
					return (1);
				else
					return (0);
			}

			//==========================================================================================

			common++;
			c++;

		} else

		{
			if (counterr >= 2) {
				flag = 1;
			}
			counterr = 0;

			for (int i = 1; i <= msearch; i++)

			{

				if ((k == 0) && (c + ora + i < s1_length)

				&& (s1[c + ora + i] == s2[c + orb + i]))

				{
					float dist = 0;
					ora = ora + i;
					orb = orb + i;
					common++;

					if (i > 1) {
						dist = i - 1 + gdist;
					} else
						dist = i + gdist;
					gdist = dist;
					if (gdist > nbdf) {

						return (0);
					}

					break;

				}

				else if ((k < 0) && (c + orb + i < s2_length)

				&& (s1[c + ora] == s2[c + orb + i]))

				{
					float dist = 0;

					if (flag == 1) {
						dist = 0;
					}

					else {
						if (i > 1) {
							dist = i - 1 + gdist;
						} else
							dist = i + gdist;
						gdist = dist;
					}
					if (gdist > nbdf) {
						return (0);
					}
					orb = orb + i;

					common++;

					break;

				}

				else if ((k > 0) && (c + ora + i < s1_length)
						&& (s1[c + ora + i] == s2[c + orb]))

						{
					float dist = 0;

					if (flag == 1) {
						dist = 0;
					}

					else {
						if (i > 1) {
							dist = i - 1 + gdist;
						} else
							dist = i + gdist;
						gdist = dist;
					}
					if (gdist > nbdf) {
						return (0);
					}

					ora = ora + i;
					common++;

					break;

				}
				if (i == msearch)
					gdist++;

			}
			c++;

		}

	}

	return (1);

}
void help() {

	printf(" \n\n ");

	printf(
			"Command line : Use parameters ( Usual cluster: From a single Input Sequence file ) : cluster max_distance , input file , output file, options \n\n");

	printf(
			"Example for Command line 1 ( Without Dereplication ) :  ./crunchclust --diff 5 --in InputSequenceFile.fasta --out OutputClusterFile.clstr --d_hl --endgaps \n\n");

	printf(
			"Example for Command line 2 ( With Dereplication ):  ./crunchclust --in InputSequenceFile.fasta --out OutputClusterFile.clstr --strict --min 10 --max 400 --kmin 0 --kmax 10 --ksteps 2 --d_hl --endgaps --\n\n");

	printf("OPTIONS for Strict Dereplication :  \n\n");

	printf("--strict   :  Whether to do strict dereplication or not\n");
	printf(
			"If you do not mention this option then crunchclust will not do any dereplication \n");
	printf(
			"After dereplication the frequency of the sequence in the raw data will be added at the end of the\n");
	printf(
			"tag name of each retained sequence which is the single representative of the\n");
	printf("duplicate  sequences. \n\n");

	printf(
			"--keep_n   :  Whether to keep the sequences with N or not during dereplication \n");
	printf(
			" If you do not mention this option then dereplication will delete the sequences containing N  \n");

	printf(
			"--min   :  Sequences shorter then this number will be deleted during dereplication, Default is 5 \n");
	printf(
			"--max   :  Sequences bigger then this number will be cropped during dereplication, Default is 400 \n");
	printf(
			" DEREPLICATION information will be printed in a file as a format called names ( it will have the name Inputefile_Dereplication.names ): mentioned by mothur website \n");

	printf(" Other options : \n\n");

	printf(
			"--diff   :  Levenstine distance threshold ( in integer number like 0, 1, 2, 3, 4 etc.)  \n\n");

	printf("--in     :  Input Sequence File  \n\n");

	printf("--out    :  Name of the outputfile  \n\n ");

	printf(
			"--kmin     :  Start clustering at minimum of this distance threshold \n\n");

	printf(
			"--kmax     :  Finish clustering at maximum of this distance threshold  \n\n");

	printf(
			"--ksteps     :  Steps of clustering at different distance threshold between --kmin to --kmax  \n\n");

	printf(
			" Two different options for distance calculations with Levinstein Distance --d_hl, --d_all \n\n ");

	printf(
			"--d_hl      :  does not count differences in ths homopolymer region, \n ");
	printf("--d_all     :  counts all the differences, \n ");

	printf(" Options for counting the end gaps \n\n ");
	printf("	              --endgaps   : For counting the end gaps   \n");
	printf(
			"	              --noendgaps : For not counting the end gaps   \n \n");
	printf(
			" PLEASE READ THE DOCUMENTATION FOR MORE DETAILS ABOUT COMMAND LINE OPTIONS \n\n ");

	exit(150);

}
//====================================================================================
int readfasta_file(ifstream &in1) {
	char c0, c1;
	int no = 0;

	c0 = '\n';
	while (1) {
		if (in1.eof())
			break;
		in1.read(&c1, 1);
		if (c1 == '>' && c0 == '\n')
			no++;
		c0 = c1;
	}
	return no;
}
//====================================================================================
int read_file_in_memory(ifstream &in1, int &SEQ_no, int *SEQ_len,
		char *SEQ_tag[], char *SEQ_final[]) {

	int MAX_SEQ = 10000;
	int MAX_DES = 1000;
	int MAX_LINE_SIZE = 10000;

	char raw_seq[MAX_SEQ], raw_des[MAX_DES];
	char buffer1[MAX_LINE_SIZE];
	raw_seq[0] = raw_des[0] = buffer1[0] = 0;
	int read_in = 0;

	int jj = -1;

	SEQ_no = 0;
	while (1) {
		if (in1.eof())
			break;
		in1.getline(buffer1, MAX_LINE_SIZE - 2, '\n');
		if (in1.fail() && (buffer1[0] == '>')) {
			printf("Problem Reading File: May be tag name is too long \n");
			help();
		}
		if (buffer1[0] == '>') {
			if (read_in) { // write previous record
				n_flag = 0;
				clean(raw_seq);
				if (strlen(raw_seq) > 5) {
					if (n_flag == 1)
						SEQ_len[SEQ_no] = 1;
					else
						SEQ_len[SEQ_no] = strlen(raw_seq);

					{

						if ((SEQ_final[SEQ_no] = new char[strlen(raw_seq) + 2])
								== NULL) {
							printf(
									"Problem Reading File: May be insufficiant memory \n");
							help();
						}

						strcpy(SEQ_final[SEQ_no], raw_seq);
					}
					SEQ_no++;
				}
			}
			strncpy(raw_des, buffer1, MAX_DES - 2);
			if ((SEQ_tag[SEQ_no] = new char[strlen(raw_des) + 2]) == NULL) {
				printf("Problem Reading File: May be insufficiant memory \n");
				help();
			}
			strcpy(SEQ_tag[SEQ_no], raw_des);

			raw_seq[0] = 0;
		} else {
			read_in = 1;
			if (strlen(raw_seq) + strlen(buffer1) >= MAX_SEQ - 1) {
				printf("Problem Reading File: May be too long sequence \n");
				help();
			}
			strcat(raw_seq, buffer1);
		}
	}

	if (1) {
		n_flag = 0;
		clean(raw_seq);

		if (strlen(raw_seq) > 5) {
			if (n_flag == 1)
				SEQ_len[SEQ_no] = 1;
			else
				SEQ_len[SEQ_no] = strlen(raw_seq);

			{

				if ((SEQ_final[SEQ_no] = new char[strlen(raw_seq) + 2]) == NULL) {
					printf("Problem Reading File: Memory problem \n");
					help();
				}

				strcpy(SEQ_final[SEQ_no], raw_seq);
			}
			SEQ_no++;
		}
	}
	in1.close();

	return 0;
}

void clean(char *seq) {
	int i, j, k;
	char c1;
	int len = strlen(seq);

	for (i = 0, j = 0; i < len; i++) {
		c1 = toupper(seq[i]);
		if (keep_n != 0) {
			if (isalpha(c1))
				seq[j++] = c1;
		} else {
			if (isalpha(c1)) {
				if (c1 == 'N') {
					n_flag = 1;
					seq[j++] = c1;
				} else
					seq[j++] = c1;
			}
		}
	}
	seq[j] = 0;
}
