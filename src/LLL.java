import java.util.ArrayList;
import static java.lang.Math.*;

public class LLL {

    private static double[]
    b1 = {1,0,0,5},
    b2 = {0,1,0,5},
    b3 = {0,0,1,5};
    private static ArrayList<double[]> b = new ArrayList<double[]>();
    private static ArrayList<double[]> b_new = new ArrayList<double[]>();

    public static void Gram_Schmidt()
    {
        double[] B = new double[b.size()];
        b_new.add(0,b.get(0));
        B[0] = pow(norm(b_new.get(0)),2);
        for(int i=1; i < 3; i++)
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

    public static void main(String[] args)
    {
        b.add(0,b1);
        b.add(1,b2);
        b.add(2,b3);

        Gram_Schmidt();

        for(int i = 0; i < b_new.size(); i++) {
            for (int j = 0; j < b_new.get(i).length; j++)
                System.out.print(b_new.get(i)[j] + " ");
            System.out.println();
        }
    }
}
