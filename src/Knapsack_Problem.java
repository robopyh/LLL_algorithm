import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.*;

public class Knapsack_Problem
{
    public static void main(String[] args)
    {
        double[] a = {12, 11, 10, 9, 33, 44, 55, 66};
        double S = 100;
        ArrayList<double[]> C = new ArrayList<>();
        ArrayList<double[]> b = new ArrayList<>();
        double[] c = new double[a.length+1];

        //step 1
        double m = ceil(sqrt(a.length)*0.5);

        //step 2
        for (int i = 0; i < a.length; i++)
        {
            c[i] = 1;
            c[a.length] = m*a[i];
            C.add(c.clone());
            Arrays.fill(c,0);
        }
        Arrays.fill(c,0.5);
        c[a.length] = S*m;
        C.add(c.clone());

        LLL lll = new LLL(C);
        b = lll.getResult();

        //step 3
    }

}
