#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>

using namespace std;

int n_nodes = 0, n_elem = 0, n_cond = 0, n_layer = 0, funcNum = 2;
double chi = 0., sigma = 0.;

struct GlobalMatrix
{
	vector<double> ggl, ggu, di;
	vector<double> l, u, d;
	vector<double> r, p, z;
	vector<double> x, temp, temp0, F;
	vector<int> ig, jg;
};

struct node
{
	double x;
	double y;
};

struct elem
{
	vector<node> knts;
	vector<int> global_num;
	int lambda, func;
};

struct condition
{
	int type;
	int vertex;
	int num;
};

struct timeMesh
{
	double t0;
	double tn;
};

double Lambda(int num, node& k)
{
	switch (num)
	{
	case 1:
		return 1;
	case 2:
		return 1;
	default:
		break;
	}
	return 1;
}

double Func(int num, node& k, elem &e, double time)
{
	switch (num)
	{
	// u(x, t) = x^2*t^3
	case(1):
		return 6 * k.x * k.x * time + 3 * k.x * k.x * time * time - 2 * time * time * time;

	// u(x, t) = x^2*t^4 гипербола
	case(2):
		return 12 * k.x * k.x * time * time + 4 * k.x * k.x * time * time * time - 2 * time * time * time * time;

	// u(x, t) = x^2*t^4 парабола
	case 3:
		return 4 * k.x * k.x * time * time * time - 2 * time * time * time * time;

	case 4:
		return sin(k.x * k.y) * (k.x * k.x + k.y * k.y + 1) + k.x * k.x * k.y - 2 * k.y;
	case 5:
		return 1;
	case 6:
		return k.x;
	case 7:
		return k.x * k.x * k.y - 2 * k.y;
	default:
		break;
	}
}

double FirstCondition(int num, node& k, double time)
{
	switch (num)
	{
	// u(x, t) = x^2*t^3
	case 1:
		return k.x * k.x * time * time * time;

	// u(x, t) = x^2*t^4 гипербола
	case 2:
		return k.x * k.x * time * time * time * time;

	// u(x, t) = x^2*t^4 парабола
	case 3:
		return k.x * k.x * time * time * time * time;;

	case 4:
		return k.x * k.x * k.y + sin(k.x * k.y);
	case 5:
		return 1;
	case 6:
		return k.x;
	case 7:
		return k.x * k.x * k.y;
	default:
		break;
	}
}

bool Input(vector<node>& nodes, vector<elem> &elems, vector<condition> &conds, timeMesh &tm)
{
	FILE* info;
	if (!fopen_s(&info, "info.txt", "r"))
	{
		fscanf_s(info, "%d %d %lf %lf", &n_nodes, &n_elem, &sigma, &chi);
	}
	else
	{
		cout << "info.txt is invalid";
		return false;
	}
	fclose(info);

	nodes.resize(n_nodes);
	elems.resize(n_elem);

	FILE* xy;
	if (!fopen_s(&xy, "xy.txt", "r"))
	{
		for (int i = 0; i < n_nodes; i++)
		{
			fscanf_s(xy, "%lf %lf", &(nodes[i].x), &(nodes[i].y));
		}
	}
	else
	{
		cout << "xy.txt is invalid";
		return false;
	}
	fclose(xy);

	FILE* el;
	if (!fopen_s(&el, "elem.txt", "r"))
	{
		for (int i = 0; i < n_elem; i++)
		{
			elems[i].global_num.resize(6);
			elems[i].knts.resize(6);
			fscanf_s(el, "%d %d %d %d %d %d ", &(elems[i].global_num[0]), &(elems[i].global_num[1]),
				&(elems[i].global_num[2]), &(elems[i].global_num[3]), &(elems[i].global_num[4]), &(elems[i].global_num[5]));
			for (int j = 0; j < 6; j++)
			{
				elems[i].global_num[j]--;
				elems[i].knts[j] = nodes[elems[i].global_num[j]];
			}
		}
	}
	else
	{
		cout << "elem.txt is invalid";
		return false;
	}
	fclose(el);

	FILE* mat;
	if (!fopen_s(&mat, "mat.txt", "r"))
	{
		for (int i = 0; i < n_elem; i++)
		{
			fscanf_s(mat, "%d %d", &(elems[i].lambda), &(elems[i].func));
		}
	}
	else
	{
		cout << "mat.txt is invalid";
		return false;
	}
	fclose(mat);

	FILE* s;
	if (!fopen_s(&s, "s1.txt", "r"))
	{
		fscanf_s(s, "%d", &n_cond);
		conds.resize(n_cond);
		for (int i = 0; i < n_cond; i++)
		{
			fscanf_s(s, "%d %d %d", &(conds[i].type), &(conds[i].num), &(conds[i].vertex));
		}
	}
	else
	{
		cout << "s1.txt is invalid";
		return false;
	}
	fclose(s);

	FILE* tMesh;
	if (!fopen_s(&tMesh, "timeMesh.txt", "r"))
	{
		fscanf_s(tMesh, "%lf %lf %d", &(tm.t0), &(tm.tn), &n_layer);
	}
	else
	{
		cout << "timeMesh.txt is invalid";
		return false;
	}
	fclose(tMesh);

	return true;
}

double DetD(elem &el)
{
	double det = (el.knts[1].x - el.knts[0].x) * (el.knts[2].y - el.knts[0].y) -
		(el.knts[2].x - el.knts[0].x) * (el.knts[1].y - el.knts[0].y);
	return det;
}

void CreateAlphas(vector<vector<double>> &alphas, elem &el, double &detD)
{
	alphas[0][0] = (el.knts[1].y - el.knts[2].y) / detD;
	alphas[0][1] = (el.knts[2].x - el.knts[1].x) / detD;
	alphas[1][0] = (el.knts[2].y - el.knts[0].y) / detD;
	alphas[1][1] = (el.knts[0].x - el.knts[2].x) / detD;
	alphas[2][0] = (el.knts[0].y - el.knts[1].y) / detD;
	alphas[2][1] = (el.knts[1].x - el.knts[0].x) / detD;
}

void addElemToGlobal(GlobalMatrix& M, int& i, int& j, double& value)
{
	int ind;

	if (i == j)
	{
		M.di[i] += value;
	}
	else
	{
		if (i > j)
		{
			for (ind = M.ig[i]; ind < M.ig[i + 1]; ind++)
			{
				if (M.jg[ind] == j)
				{
					M.ggl[ind] += value;
					break;
				}
			}
		}
		else
		{
			for (ind = M.ig[j]; ind < M.ig[j + 1]; ind++)
			{
				if (M.jg[ind] == i)
				{
					M.ggu[ind] += value;
					break;
				}
			}
		}
	}
}

//void CreateGlobalSystem(vector<vector<double>>& A, vector<double>& b, vector<localMatrix>& arrA, vector<vector<double>>& arrB, vector<elem>& elems)
//{
//	for (int k = 0; k < n_elem; k++)
//	{
//		for (int i = 0; i < 6; i++)
//		{
//			int indexB = elems[k].global_num[i];
//			int indexI = elems[arrA[k].elemNum].global_num[i];
//			for (int j = 0; j < 6; j++)
//			{
//				int indexJ = elems[arrA[k].elemNum].global_num[j];
//				A[indexI][indexJ] += arrA[k].localM[i][j];
//			}
//			b[indexB] += arrB[k][i];
//		}
//	}
//}

void addLocalToGlobal(GlobalMatrix& M, elem& el, vector<vector<double>>& localM)
{
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			addElemToGlobal(M, el.global_num[i], el.global_num[j], localM[i][j]);
		}
	}
}

void addLocalBToGlobal(GlobalMatrix& M, elem& el, vector<double>& localB)
{
	for (int i = 0; i < 6; i++)
	{
		M.F[el.global_num[i]] += localB[i];
	}
}

void CalcLocal(GlobalMatrix &m, GlobalMatrix& g, elem &el)
{
	vector<vector<double>> M(6), G(6), alphas(3);
	double detD = 0, buff = 0., l1 = 0., l2 = 0., l3 = 0.;
	for (int i = 0; i < 6; i++)
	{
		M[i].resize(6);
		G[i].resize(6);
	}

	for (int i = 0; i < 3; i++)
	{
		alphas[i].resize(2);
	}

	detD = DetD(el);
	CreateAlphas(alphas, el, detD);

	detD = abs(detD);
	buff = 1;
	M[0][0] = (1. / 60) * buff;
	M[1][1] = (1. / 60) * buff;
	M[2][2] = (1. / 60) * buff;
	M[3][3] = (4. / 45) * buff;
	M[4][4] = (4. / 45) * buff;
	M[5][5] = (4. / 45) * buff;
	M[0][1] = M[1][0] = (-1. / 360) * buff;
	M[0][2] = M[2][0] = (-1. / 360) * buff;
	M[1][2] = M[2][1] = (-1. / 360) * buff;
	M[3][4] = M[4][3] = (2. / 45) * buff;
	M[3][5] = M[5][3] = (2. / 45) * buff;
	M[4][5] = M[5][4] = (2. / 45) * buff;
	M[2][3] = M[3][2] = M[0][4] = M[4][0] = M[1][5] = M[5][1] = (-1. / 90) * buff;

	l1 = Lambda(el.lambda, el.knts[0]);
	l2 = Lambda(el.lambda, el.knts[1]);
	l3 = Lambda(el.lambda, el.knts[2]);

	G[0][0] = (alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (3 * l1 + l2 + l3) / 10;
	G[1][1] = (alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (l1 + 3 * l2 + l3) / 10;
	G[2][2] = (alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (l1 + l2 + 3 * l3) / 10;

	G[3][3] = (4. / 15) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (l1 + 3 * l2 + l3) +
		(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (2 * l1 + 2 * l2 + l3) +
		(alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (3 * l1 + l2 + l3));

	G[4][4] = (4. / 15) * ((alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (l1 + l2 + 3 * l3) +
		(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + 2 * l2 + 2 * l3) +
		(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (l1 + 3 * l2 + l3));

	G[5][5] = (4. / 15) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (l1 + l2 + 3 * l3) +
		(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (2 * l1 + l2 + 2 * l3) +
		(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (3 * l1 + l2 + l3));

	G[0][1] = G[1][0] = -(1. / 30) * (alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (2 * l1 + 2 * l2 + l3);
	G[0][2] = G[2][0] = -(1. / 30) * (alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (2 * l1 + l2 + 2 * l3);
	G[1][2] = G[2][1] = -(1. / 30) * (alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + 2 * l2 + 2 * l3);

	G[3][4] = G[4][3] = (4. / 15) * ((alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (l1 / 2 + l2 + l3) +
		(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (l1 + 3 * l2 + l3) +
		(alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (l1 + l2 / 2 + l3) +
		(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + l2 + l3 / 2));

	G[3][5] = G[5][3] = (4. / 15) * ((alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (l1 + l2 + l3 / 2) +
		(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (l1 + l2 / 2 + l3) +
		(alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (l1 / 2 + l2 + l3) +
		(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (3 * l1 + l2 + l3));

	G[4][5] = G[5][4] = (4. / 15) * ((alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (l1 / 2 + l2 + l3) +
		(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (l1 + l2 / 2 + l3) +
		(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (l1 + l2 + l3 / 2) +
		(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (l1 + l2 + 3 * l3));

	G[0][3] = G[3][0] = (1. / 30) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (3 * l1 - 2 * l2 - l3) +
		(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (14 * l1 + 3 * l2 + 3 * l3));

	G[1][3] = G[3][1] = (1. / 30) * ((alphas[1][0] * alphas[1][0] + alphas[1][1] * alphas[1][1]) * (-2 * l1 + 3 * l2 - l3) +
		(alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (3 * l1 + 14 * l2 + 3 * l3));

	G[2][3] = G[3][2] = (1. / 30) * ((alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (-2 * l1 - l2 + 3 * l3) +
		(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (-1 * l1 - 2 * l2 + 3 * l3));

	G[0][4] = G[4][0] = (1. / 30) * ((alphas[1][0] * alphas[0][0] + alphas[0][1] * alphas[1][1]) * (3 * l1 - l2 - 2 * l3) +
		(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (3 * l1 - 2 * l2 - l3));

	G[1][4] = G[4][1] = (1. / 30) * ((alphas[1][0] * alphas[0][0] + alphas[1][1] * alphas[0][1]) * (-1 * l1 + 3 * l2 - 2 * l3) +
		(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (3 * l1 + 14 * l2 + 3 * l3));

	G[2][4] = G[4][2] = (1. / 30) * ((alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (3 * l1 + 3 * l2 + 14 * l3) +
		(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (-1 * l1 - 2 * l2 + 3 * l3));

	G[0][5] = G[5][0] = (1. / 30) * ((alphas[0][0] * alphas[0][0] + alphas[0][1] * alphas[0][1]) * (3 * l1 - l2 - 2 * l3) +
		(alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (14 * l1 + 3 * l2 + 3 * l3));

	G[1][5] = G[5][1] = (1. / 30) * ((alphas[0][0] * alphas[1][0] + alphas[0][1] * alphas[1][1]) * (-1 * l1 + 3 * l2 - 2 * l3) +
		(alphas[1][0] * alphas[2][0] + alphas[1][1] * alphas[2][1]) * (-2 * l1 + 3 * l2 - l3));

	G[2][5] = G[5][2] = (1. / 30) * ((alphas[0][0] * alphas[2][0] + alphas[0][1] * alphas[2][1]) * (3 * l1 + 3 * l2 + 14 * l3) +
		(alphas[2][0] * alphas[2][0] + alphas[2][1] * alphas[2][1]) * (-2 * l1 - l2 + 3 * l3));

	for (int k = 0; k < 6; k++)
	{
		for (int z = 0; z < 6; z++)
		{
			M[k][z] = detD * M[k][z];
			G[k][z] = detD * G[k][z];
		}
	}

	addLocalToGlobal(m, el, M);
	addLocalToGlobal(g, el, G);
}

void CalcLocalB(GlobalMatrix &A, elem &el, double time)
{
	double buff = 0., f1 = 0., f2 = 0., f3 = 0., f4 = 0., f5 = 0., f6 = 0.;
	vector<double> arr(6);

	double detD = abs(DetD(el));
	buff = detD / 360;

	f1 = Func(el.func, el.knts[0], el, time);
	f2 = Func(el.func, el.knts[1], el, time);
	f3 = Func(el.func, el.knts[2], el, time);
	f4 = Func(el.func, el.knts[3], el, time);
	f5 = Func(el.func, el.knts[4], el, time);
	f6 = Func(el.func, el.knts[5], el, time);

	arr[0] = buff * (6 * f1 - f2 - f3 - 4 * f5);
	arr[1] = buff * (-f1 + f2 * 6 - f3 - 4 * f6);
	arr[2] = buff * (-f1 - f2 + 6 * f3 - 4 * f4);

	buff = detD / 45;

	arr[3] = buff * (-f3 / 2 + 4 * f4 + 2 * f5 + 2 * f6);
	arr[4] = buff * (-f1 / 2 + 2 * f4 + 4 * f5 + 2 * f6);
	arr[5] = buff * (-f2 / 2 + 2 * f4 + 2 * f5 + 4 * f6);

	addLocalBToGlobal(A, el, arr);
}

//void AccountConditions_1(vector<vector<double>> &A, vector<double> &b, vector<condition> &conds, vector<node> &nodes)
//{
//	for (int i = 0; i < n_cond; i++)
//	{
//		int index = conds[i].vertex - 1;
//		for (int j = 0; j < n_nodes; j++)
//		{
//			A[index][j] = (index == j ? 1 : 0);
//		}
//		b[index] = FirstCondition(conds[i].num, nodes[index]);
//	}
//}

void AccountConditions(GlobalMatrix& A, vector<double>& b, vector<condition>& conds, vector<node>& nodes, double &time)
{
	int ibeg, iend, jind;
	for (int i = 0; i < n_cond; i++)
	{
		int index = conds[i].vertex - 1;
		A.di[index] = 1;
		b[index] = FirstCondition(conds[i].num, nodes[index], time);

		ibeg = A.ig[index];
		iend = A.ig[index + 1];
		for (int j = ibeg; j < iend; j++)
		{
			A.ggl[j] = 0;
		}
		// Обнуляем внедиагональные элементы i-ой строки в ggu
		jind = index;
		for (int j = 0; j < A.ig[n_nodes]; j++)
		{
			if (A.jg[j] == jind)
			{
				A.ggu[j] = 0;
			}
		}
	}
}

//double ScalarMult(vector<double> &vec1, vector<double>& vec2)
//{
//	double res = 0.;
//	for (int i = 0; i < n_nodes; i++)
//	{
//		res += vec1[i] * vec2[i];
//	}
//	return res;
//}
//
//void MatrixVectorMult(vector<vector<double>>& A, vector<double> &vec, vector<double> &res)
//{
//	for (int i = 0; i < n_nodes; i++)
//	{
//		res[i] = ScalarMult(A[i], vec);
//	}
//}
//
//void LOC(vector<vector<double>> &A, vector<double> &b, vector<double> &x)
//{
//	double residual = 0., residual1 = 0., alpha = 0., beta = 0., scMult = 0., eps = 1e-15;
//	vector<double> r(n_nodes), z(n_nodes), p(n_nodes);
//
//	for (int i = 0; i < n_nodes; i++)
//	{
//		b[i] -= ScalarMult(A[i], x);
//	}
//	r = z = b;
//	MatrixVectorMult(A, z, p);
//	residual = ScalarMult(r, r);
//
//	for (int i = 1; i < 100000 && residual > eps && residual != residual1; i++)
//	{
//		residual1 = residual;
//		scMult = ScalarMult(p, p);
//		alpha = ScalarMult(p, r) / scMult;
//
//		for (int j = 0; j < n_nodes; j++)
//		{
//			x[j] += alpha * z[j];
//			r[j] -= alpha * p[j];
//		}
//
//		MatrixVectorMult(A, r, b);
//
//		beta = -ScalarMult(p, b) / scMult;
//
//		for (int j = 0; j < n_nodes; j++)
//		{
//			z[j] = r[j] + beta * z[j];
//			p[j] = b[j] + beta * p[j];
//		}
//
//		residual -= alpha * alpha * scMult;
//	}
//}

// Процедура построения портрета матрицы
void Portrait(GlobalMatrix &M, vector<elem> &elems)
{
	vector<int> listbeg(n_nodes);
	vector<int> list1(8 * n_nodes);
	vector<int> list2(8 * n_nodes);
	
	int ielem, i, k, j;
	int ind1, ind2, iaddr;
	int listsize = -1;

	for (ielem = 0; ielem < n_elem; ielem++)
	{
		for (i = 0; i < 6; i++)
		{
			k = elems[ielem].global_num[i];
			for (j = i + 1; j < 6; j++)
			{
				ind1 = k;
				ind2 = elems[ielem].global_num[j];
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				iaddr = listbeg[ind2];
				if (iaddr == 0)
				{
					listsize++;
					listbeg[ind2] = listsize;
					list1[listsize] = ind1;
					list2[listsize] = 0;
				}
				else
				{
					while (list1[iaddr] < ind1 && list2[iaddr] > 0)
					{
						iaddr = list2[iaddr];
					}
					if (list1[iaddr] > ind1)
					{
						listsize++;
						list1[listsize] = list1[iaddr];
						list2[listsize] = list2[iaddr];
						list1[iaddr] = ind1;
						list2[iaddr] = listsize;
					}
					else
					{
						if (list1[iaddr] < ind1)
						{
							listsize++;
							list2[iaddr] = listsize;
							list1[listsize] = ind1;
							list2[listsize] = 0;
						}
					}
				}
			}
		}
	}

	M.ig.resize(n_nodes + 1);
	M.di.resize(n_nodes);
	M.ggl.resize(list1.size());
	M.ggu.resize(list1.size());
	M.u.resize(list1.size());
	M.l.resize(list1.size());
	M.d.resize(n_nodes);
	M.jg.resize(list1.size());
	M.x.resize(n_nodes);
	M.temp.resize(n_nodes);
	M.temp0.resize(n_nodes);
	M.F.resize(n_nodes);
	M.r.resize(n_nodes);
	M.p.resize(n_nodes);
	M.z.resize(n_nodes);

	M.ig[0] = 0;

	for (i = 0; i < n_nodes; i++)
	{
		M.ig[i + 1] = M.ig[i];
		iaddr = listbeg[i];
		while (iaddr != 0)
		{
			M.jg[M.ig[i + 1] + 1] = list1[iaddr];
			M.ig[i + 1]++;
			iaddr = list2[iaddr];
		}
	}
	for (i = 3; i < n_nodes + 1; i++)
	{
		M.ig[i]++;
	}
}

void multyMatrixVector(GlobalMatrix &m, vector<double> &x, vector<double> &res)
{
	for (int i = 0; i < n_nodes; ++i) {
		int gi = m.ig[i], gi_1 = m.ig[i + 1];
		res[i] = m.di[i] * x[i];
		for (int j = gi; j < gi_1; ++j) {
			int column = m.jg[j];
			res[i] += m.ggl[j] * x[column];
			res[column] += m.ggu[j] * x[i];
		}
	}
}

double scal(vector<double> &a, vector<double> &b)
{
	double s = 0.0;
	for (int i = 0; i < n_nodes; i++)
		s += a[i] * b[i];
	return s;
}

double norm(vector<double> &a)
{
	return sqrt(scal(a, a));
}

void LOS(GlobalMatrix &m)
{
	int count = 0, maxiter = 100000;
	double eps = 1e-20;

	for (int i = 0; i < n_nodes; ++i)
	{
		m.x[i] = 1;
	}

	multyMatrixVector(m, m.x, m.temp);

	for (int i = 0; i < n_nodes; ++i)
	{
		m.r[i] = m.F[i] - m.temp[i];
		m.z[i] = m.r[i];
	}

	multyMatrixVector(m, m.z, m.p);
	double sr = scal(m.r, m.r);

	while (sr > eps && count <= maxiter)
	{
		double pp = scal(m.p, m.p);

		double ak = scal(m.p, m.r) / pp;

		for (int i = 0; i < n_nodes; ++i)
		{
			m.x[i] = m.x[i] + ak * m.z[i];
			m.r[i] = m.r[i] - ak * m.p[i];
		}

		multyMatrixVector(m, m.r, m.temp);

		double bk = -scal(m.p, m.temp) / pp;

		for (int i = 0; i < n_nodes; ++i)
		{
			m.z[i] = m.r[i] + bk * m.z[i];
			m.p[i] = m.temp[i] + bk * m.p[i];
		}
		sr = sqrt(scal(m.r, m.r));
		++count;
	}
}

void rndPoint(node &k, elem &el, vector<double> &q, double &time)
{
	vector<vector<double>> alpha(3);
	for (int i = 0; i < 3; i++)
		alpha[i].resize(3);
	vector<double> L(3);
	vector<double> basis(6);

	double detD = DetD(el);
	double S23 = abs((el.knts[2].x - el.knts[1].x) * (k.y - el.knts[1].y) - (k.x - el.knts[1].x) * (el.knts[2].y - el.knts[1].y));
	double S31 = abs((el.knts[0].x - el.knts[2].x) * (k.y - el.knts[2].y) - (k.x - el.knts[2].x) * (el.knts[0].y - el.knts[2].y));
	double S12 = abs((el.knts[1].x - el.knts[0].x) * (k.y - el.knts[0].y) - (k.x - el.knts[0].x) * (el.knts[1].y - el.knts[0].y));

	L[0] = S23 / abs(detD);
	L[1] = S31 / abs(detD);
	L[2] = S12 / abs(detD);

	for (int i = 0; i < 3; i++)
		basis[i] = L[i] * (2. * L[i] - 1);
	basis[3] = 4. * L[0] * L[1];
	basis[4] = 4. * L[1] * L[2];
	basis[5] = 4. * L[0] * L[2];
	double sum = 0;
	// для вершин
	sum += q[el.global_num[0]] * basis[0];
	sum += q[el.global_num[1]] * basis[1];
	sum += q[el.global_num[2]] * basis[2];

	// для середин
	sum += q[el.global_num[3]] * basis[3];
	sum += q[el.global_num[4]] * basis[4];
	sum += q[el.global_num[5]] * basis[5];
	double resh = FirstCondition(3, k, time);
	cout << setprecision(3) << sum << "\t" << "\t" << resh << "\t" << "\t" << fabs(resh - sum) << "\n\n";
}

void calcGlobalMandG(vector<elem>& elems, GlobalMatrix &M, GlobalMatrix &G)
{
	for (int i = 0; i < n_elem; i++)
	{
		CalcLocal(M, G, elems[i]);
	}
}

void calcB(vector<elem>& elems, GlobalMatrix& A, double &time)
{
	for (int i = 0; i < n_elem; i++)
	{
		CalcLocalB(A, elems[i], time);
	}
}

void calcFirstThreeQs(timeMesh &tm, vector<vector<double>> &solution, vector<node> &k)
{
	double h = (tm.tn - tm.t0) / n_layer, currT = tm.t0;
	for (int i = 0; i < 3; i++)
	{
		currT = tm.t0 + i * h;
		for (int j = 0; j < n_nodes; j++)
			solution[i][j] = FirstCondition(funcNum, k[j], currT);
	}
}

void writeResults(string fname, vector<vector<double>>& solutions, vector<node> &nodes, timeMesh &tm)
{
	ofstream file(fname);
	double h = (tm.tn - tm.t0) / n_layer, currT = tm.t0;

	for (int i = 0; i < solutions.size(); i++)
	{
		currT = tm.t0 + i * h;
		file << "q" << i << endl;
		file << "calc\t\t\t" << "analitic\t\t" << "diff" << endl;

		for (int j = 0; j < n_nodes; j++)
		{
			double a = FirstCondition(funcNum, nodes[j], currT);
			file << scientific << solutions[i][j] << "\t" << a << "\t" << abs(solutions[i][j] - a) << endl;
		}
		file << endl << endl;
	}
}

int main()
{
	GlobalMatrix M;
	vector<node> nodes;
	vector<elem> elems;
	vector<condition> conds;
	timeMesh tm;

	if (!Input(nodes, elems, conds, tm)) return 1;
	double h = (tm.tn - tm.t0) / n_layer, currT = tm.t0;

	vector<vector<double>> solutions(n_layer + 1);
	for (int i = 0; i <= n_layer; i++)
	{
		solutions[i].resize(n_nodes);
	}

	Portrait(M, elems);
	GlobalMatrix G = M, A = M, T = M;

	calcGlobalMandG(elems, M, G);
	calcFirstThreeQs(tm, solutions, nodes);
 
	//vector<vector<double>> A_1(n_nodes);
	//vector<double> b(n_nodes);
	/*for (int i = 0; i < n_nodes; i++)
		A_1[i].resize(n_nodes);*/

	//CreateGlobalSystem(A_1, b, arrayLocals, arrayLocalsB, elems);
	double t01 = 0, t02 = 0, t03 = 0, t10 = 0, t12 = 0, t13 = 0, t20 = 0, t21 = 0,
		t23 = 0, t30 = 0, t31 = 0, t32 = 0, t_chi_j = 0, t_sigma_j = 0, t_chi_j3 = 0,
		t_chi_j2 = 0, t_chi_j1 = 0, t_sigma_j3 = 0, t_sigma_j2 = 0, t_sigma_j1 = 0;

	vector<double> temp(n_nodes);

	for (int i = 3; i <= n_layer; i++)
	{
		// Текущий временной слой
		currT = tm.t0 + i * h;
		// Изменения по t
		t01 = currT - (currT - h);
		t02 = currT - (currT - 2 * h);
		t03 = currT - (currT - 3 * h);
		t10 = (currT - h) - currT;
		t12 = (currT - h) - (currT - 2 * h);
		t13 = (currT - h) - (currT - 3 * h);
		t20 = (currT - 2 * h) - currT;
		t21 = (currT - 2 * h) - (currT - h);
		t23 = (currT - 2 * h) - (currT - 3 * h);
		t30 = (currT - 3 * h) - currT;
		t31 = (currT - 3 * h) - (currT - h);
		t32 = (currT - 3 * h) - (currT - 2 * h);

		// Вектор правой части на текущем временном слое
		A.F.assign(A.F.size(), 0);
		calcB(elems, A, currT);

		// Множители для неявной 4-х слойной схемы
		t_chi_j = 2 * (t01 + t02 + t03) / (t01 * t02 * t03) * chi;
		t_sigma_j = (t02 * t01 + t03 * t02 + t01 * t03) / (t03 * t02 * t01) * sigma;
		t_chi_j3 = 2 * (t01 + t02) / (t32 * t31 * t30) * chi;
		t_sigma_j3 = (t02 * t01) / (t32 * t31 * t30) * sigma;
		t_chi_j2 = 2 * (t01 + t03) / (t23 * t21 * t20) * chi;
		t_sigma_j2 = (t01 * t03) / (t23 * t21 * t20) * sigma;
		t_chi_j1 = 2 * (t02 + t03) / (t13 * t12 * t10) * chi;
		t_sigma_j1 = (t03 * t02) / (t13 * t12 * t10) * sigma;

		// Левая часть
		for (int j = 0; j < A.ig[n_nodes - 1]; j++)
		{
			A.ggl[j] = t_chi_j * M.ggl[j] + t_sigma_j * M.ggl[j] + G.ggl[j];
			A.ggu[j] = t_chi_j * M.ggu[j] + t_sigma_j * M.ggu[j] + G.ggu[j];
		}

		for (int j = 0; j < n_nodes; j++)
		{
			A.di[j] = t_chi_j * M.di[j] + t_sigma_j * M.di[j] + G.di[j];
		}

		// Правая часть
		// deltaT * M * q[j-3]
		for (int j = 0; j < A.ig[n_nodes - 1]; j++)
		{
			T.ggl[j] = t_chi_j3 * M.ggl[j] + t_sigma_j3 * M.ggl[j];
			T.ggu[j] = t_chi_j3 * M.ggu[j] + t_sigma_j3 * M.ggu[j];
		}

		for (int j = 0; j < n_nodes; j++)
		{
			T.di[j] = t_chi_j3 * M.di[j] + t_sigma_j3 * M.di[j];
		}

		multyMatrixVector(T, solutions[i - 3], temp);

		for (int j = 0; j < n_nodes; j++)
		{
			A.F[j] -= temp[j];
		}

		// deltaT * M * q[j-2]
		for (int j = 0; j < A.ig[n_nodes - 1]; j++)
		{
			T.ggl[j] = t_chi_j2 * M.ggl[j] + t_sigma_j2 * M.ggl[j];
			T.ggu[j] = t_chi_j2 * M.ggu[j] + t_sigma_j2 * M.ggu[j];
		}

		for (int j = 0; j < n_nodes; j++)
		{
			T.di[j] = t_chi_j2 * M.di[j] + t_sigma_j2 * M.di[j];
		}

		multyMatrixVector(T, solutions[i - 2], temp);

		for (int j = 0; j < n_nodes; j++)
		{
			A.F[j] -= temp[j];
		}

		// deltaT * M * q[j-1]
		for (int j = 0; j < A.ig[n_nodes - 1]; j++)
		{
			T.ggl[j] = t_chi_j1 * M.ggl[j] + t_sigma_j1 * M.ggl[j];
			T.ggu[j] = t_chi_j1 * M.ggu[j] + t_sigma_j1 * M.ggu[j];
		}

		for (int j = 0; j < n_nodes; j++)
		{
			T.di[j] = t_chi_j1 * M.di[j] + t_sigma_j1 * M.di[j];
		}

		multyMatrixVector(T, solutions[i - 1], temp);

		for (int j = 0; j < n_nodes; j++)
		{
			A.F[j] -= temp[j];
		}

		AccountConditions(A, A.F, conds, nodes, currT);
		LOS(A);
		solutions[i].swap(A.x);
	}

	string fname = to_string(funcNum) + "_" + to_string(n_layer) + ".txt";

	writeResults(fname, solutions, nodes, tm);

	//AccountConditions_1(A_1, b, conds, nodes);
	//LOC(A_1, b, x);

	/*node k;
	k.x = 1;
	k.y = 0.25;*/

	//rndPoint(k, elems[0], A.x);

	
	
	return 0;
}