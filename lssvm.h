#pragma once
#include "common.h"


class lssvm
{

public:
    long feature_dim;
    long sample_size;
    vector<vector<double>> x_train;
    vector<double> y_train;
    double learning_rate = 0.01;
    double constraint = 0.1;
    double gamma = 0.1;
    vector<double> weight;
    vector<double> tilda_one;
    vector<vector<double>> linearMat; // (sample_size+1) x (sample_size+1)
    string tag = "linear";

    lssvm(vector<vector<double>> &train_x, vector<double> &train_y, long featrue_dim, long sample_size,double constraint = 1.0,double gamma = 0.1)
    {

        this->constraint = constraint;

        x_train.resize(sample_size);
        for (int i = 0; i < sample_size; ++i)
        {
            x_train[i].resize(featrue_dim);
        }

        y_train.resize(sample_size);

        weight.resize(sample_size + 1);
        for (int i = 0; i <= sample_size; ++i)
        {
            weight[i] = 0.1;
        }

        tilda_one.resize(sample_size + 1);
        tilda_one[0] = 0.0;

        for (int i = 1; i < sample_size + 1; ++i)
        {
            tilda_one[i] = 1.0;
        }

        x_train = train_x;
        // cout << 2;
        y_train = train_y;
        // cout << 3;
        this->feature_dim = featrue_dim;
        // cout << 4;
        this->sample_size = sample_size;
        // cout << 5;

        gen_linear(gamma);

        fit();

        validation();
    }

    vector<double> predict(vector<vector<double>> &test_x, long test_size)
    {   


        vector<double> rtn(test_size);

        for (int i = 0; i < test_size; ++i)
        {
            double prob = 0.0;

            for (int j = 0; j < sample_size; ++j)
            {
                vector<double> tmp = test_x[i];

                // cout << "rbf\n";
                double kkk = 0.0;

                //sub(tmp, x_train[j], tmp);
                for (int k = 0; k < feature_dim; ++k)
                {
                    tmp [k]= tmp[k] - x_train[j][k];
                }
                
                for (int k = 0; k < feature_dim; ++k)
                {
                    kkk += tmp[k] * tmp[k];
                }

                // cout << "prob: " << prob << '\n';

                kkk = exp(-kkk);

                prob += weight[j+1] * y_train[j] * kkk;
            }

            prob += weight[0];

            // cout << "signal: " << prob << '\n';

            prob = 1.0 / (1 + exp(-prob));

            // cout << "prob: " << prob << '\n';

            rtn[i] = (prob > 0.5) ? 1.0 : -1.0;
        }
        
        return rtn;
    }

    void show_linear()
    {

        cout << "Linear matrix:\n";
        for (int i = 0; i < sample_size + 1; ++i)
        {
            for (int j = 0; j < sample_size + 1; ++j)
            {
                cout << linearMat[i][j] << " ";
            }
            cout << '\n';
        }
        cout << '\n';
    }

    void show_weight()
    {
        cout << "initial weight: \n";
        for (int j = 0; j < sample_size + 1; ++j)
        {
            cout << weight[j] << " ";
        }
        cout << '\n';
    }

    void show_tilda()
    {
        cout << "tilda one!: \n";
        for (int j = 0; j < sample_size + 1; ++j)
        {
            cout << tilda_one[j] << " ";
        }
        cout << '\n';
    }

    void gen_kernel(double gamma)
    {
        for (int i = 0; i < sample_size; ++i)
        {
            for (int j = 0; j < sample_size; ++j)
            {
                double tmp = 0.0;
                for (int k = 0; k < feature_dim; ++k)
                {
                    tmp += (x_train[i][k] - x_train[j][k]) * (x_train[i][k] - x_train[j][k]);
                }
                linearMat[i + 1][j + 1] = exp(-gamma*tmp);
                linearMat[i + 1][j + 1] = y_train[i] * y_train[j] * linearMat[i + 1][j + 1];
            }
        }
        for (int i = 0; i < sample_size; ++i)
        {
            linearMat[i + 1][i + 1] += 0.1;
        }
    }

    void gen_linear(double gamma)
    {

        linearMat.resize(sample_size + 1);
        for (int i = 0; i < sample_size + 1; ++i)
        {
            linearMat[i].resize(sample_size + 1);
        }

        gen_kernel(gamma);
        linearMat[0][0] = 0.0;

        for (int i = 1; i <= sample_size; ++i)
        {
            linearMat[0][i] = -y_train[i - 1];
            linearMat[i][0] = y_train[i - 1];
        }
    }

    void fit()
    {

        double current_loss = 100.0;
        for (int i = 1; i <= 100000; ++i)
        {
            vector<double> tmp = weight;

            mult(linearMat, weight, weight);
            sub(weight, tilda_one, weight);
            mult(linearMat, weight, weight);
            constmult(weight, learning_rate, weight);
            sub(tmp, weight, weight);

            if (i % 1000 == 0)
            {
                cout << i << "-th repeat-> ";
                validation();
                cout << '\n';
            }
        }
    }

    void validation()
    {
        vector<double> tmp = weight;

        mult(linearMat, weight, tmp);

        sub(tmp, tilda_one, tmp);

        double var = 0.0;
        for (int i = 0; i < sample_size + 1; ++i)
        {
            // cout << abs(tmp[i]) << " ";
            var += tmp[i] * tmp[i];
        }

        cout << "norm: " << var / double(sample_size + 1) << '\n';
        // cout << "loss: " << var/(sample_size+1) << '\n';
    }

    void mult(vector<vector<double>> &mat, vector<double> &vec, vector<double> &rtn)
    {
        vector<double> t = rtn;
        for (int i = 0; i < sample_size + 1; ++i)
        {
            double tmp = 0.0;
            for (int j = 0; j < sample_size + 1; ++j)
            {
                tmp += (mat[i][j] * vec[j]);
            }
            t[i] = tmp;
        }
        rtn = t;
    }

    void add(vector<double> &vec1, vector<double> &vec2, vector<double> &rtn)
    {
        vector<double> t = rtn;
        for (int i = 0; i < sample_size + 1; ++i)
        {
            t[i] = vec1[i] + vec2[i];
        }
        rtn = t;
    }

    void sub(vector<double> &vec1, vector<double> &vec2, vector<double> &rtn)
    {
        vector<double> t = rtn;
        for (int i = 0; i < sample_size + 1; ++i)
        {
            t[i] = vec1[i] - vec2[i];
        }
        rtn = t;
    }

    void constmult(vector<double> &vec, double alpha, vector<double> &rtn)
    {
        vector<double> t = rtn;
        for (int i = 0; i < sample_size + 1; ++i)
        {
            t[i] = alpha * vec[i];
        }
        rtn = t;
    }
};