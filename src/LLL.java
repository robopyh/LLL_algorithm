import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import static java.lang.Math.*;

public class LLL {

    private static ArrayList<double[]> b;
    private static ArrayList<double[]> b_new = new ArrayList<>();
    private static double[] B;

    private static void Gram_Schmidt()
    {
        B = new double[b.size()];
        b_new.add(0,b.get(0));
        B[0] = pow(norm(b_new.get(0)),2);
        for(int i=1; i < b.size(); i++)
        {
            b_new.add(i,b.get(i));
            for(int j=0; j < i; j++)
            {
                double mu = scalar_multiply(b.get(i), b_new.get(j)) / B[j];
                double[] b_temp = subtraction(b_new.get(i), scalar_multiply(mu, b_new.get(j)));
                b_new.set(i,b_temp);
            }
            B[i] = pow(norm(b_new.get(i)),2);
        }
    }

    private static double norm(double[] vector)
    {
        return sqrt(scalar_multiply(vector, vector));
    }

    private static double[] scalar_multiply(double c, double[] b)
    {
        double[] result = new double[b.length];
        for(int i=0; i < b.length; i++)
            result[i] = c * b[i];
        return result;
    }

    private static double scalar_multiply(double[] a, double[] b)
    {
        double result = 0;
        for(int i=0; i < a.length; i++)
            result += a[i] * b[i];
        return result;
    }

    private static double[] subtraction(double[] a, double[] b)
    {
        double[] result = new double[a.length];
        for(int i=0; i < a.length; i++)
            result[i] = a[i] - b[i];
        return result;
    }

    private static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    private static void LLL_algorithm()
    {
        //step 1
        Gram_Schmidt();

        //debug
        //step 2
        int k = 1;
        double r;
        double mu_temp;
        double[][] mu = new double[b.size()][b.size()];

        for (int i = 0; i < b.size(); i++)
            for (int j = 0; j < b.size(); j++)
                mu[i][j] = scalar_multiply(b.get(i), b_new.get(j)) / B[j];

        while(true)
        {
            //step 3
            if (abs(mu[k][k - 1]) > 0.5)
            {
                //step 3.1
                r = (mu[k][k - 1] > 0) ? floor(0.5 + mu[k][k - 1]) : -floor(0.5 - mu[k][k - 1]);
                //step 3.2
                b.set(k, subtraction(b.get(k), scalar_multiply(r, b.get(k - 1))));
                //step 3.3
                for (int j = 0; j < k - 1; j++)
                    mu[k][j] = mu[k][j] - r * mu[k - 1][j];
                //step 3.4
                mu[k][k - 1] = mu[k][k - 1] - r;
            }

            //step 4
            if (B[k] < (0.75 - pow(mu[k][k - 1], 2)) * B[k - 1]) {
                //step 4.1
                mu_temp = mu[k][k - 1];
                double B_temp = B[k] + pow(mu_temp, 2) * B[k - 1];
                mu[k][k - 1] = mu_temp * B[k - 1] / B_temp;
                B[k] = B[k - 1] * B[k] / B_temp;
                B[k - 1] = B_temp;
                //step 4.2
                double[] temp = b.get(k);
                b.set(k, b.get(k - 1));
                b.set(k - 1, temp);
                //step 4.3
                if (k > 1)
                    for (int j = 0; j < k - 1; j++) {
                        mu_temp = mu[k][j];
                        mu[k][j] = mu[k - 1][j];
                        mu[k - 1][j] = mu_temp;
                    }
                //step 4.4
                for (int s = k + 1; s < b.size(); s++) {
                    double t = mu[s][k];
                    mu[s][k] = mu[s][k - 1] - mu_temp * t;
                    mu[s][k - 1] = t + mu[k][k - 1] * mu[s][k];
                }
                //step 4.5
                k = max(1, k - 1);
            }
            //step 5
            else
            {
                //step 5.1
                for (int l = k-1; l >= 0; l--)
                    if(abs(mu[k][l]) > 0.5)
                    {
                        //step 5.1.1
                        r = (mu[k][l] > 0) ? floor(0.5 + mu[k][l]) : floor(-(0.5 + mu[k][l]));
                        //step 5.1.2
                        b.set(k, subtraction(b.get(k), scalar_multiply(r, b.get(l))));
                        //step 5.1.3
                        for(int j = 0; j < l; j++)
                            mu[k][j] = mu[k][j] - r*mu[l][j];
                        //step 5.1.4
                        mu[k][l] = mu[k][l] - r;
                    }
                //step 5.2
                    k++;
            }

            if(k > b.size() - 1)
                break;
        }

        //LLL result
        for (double[] aB : b) {
            for (double anAB : aB) System.out.print(anAB + " ");
            System.out.println();
        }
    }

    ArrayList<double[]> getResult()
    {
        return b;
    }

    public LLL(ArrayList<double[]> vectors)
    {
        /*double[] b1 = {1,0,1,2};
        double[] b2 = {1,-1,2,0};
        double[] b3 = {-1,2,0,1};
        b = new ArrayList<>();
        b.add(b1);
        b.add(b2);
        b.add(b3);*/
        b = new ArrayList<>(vectors);
        LLL_algorithm();
    }
}
