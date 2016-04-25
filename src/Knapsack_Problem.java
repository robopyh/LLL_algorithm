import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.*;

public class Knapsack_Problem
{
    private static BigDecimal zero = new BigDecimal("0");
    private static BigDecimal half = new BigDecimal("0.5");

    public static void main(String[] args)
    {
        int[] a = {12, 11, 10, 9, 33, 44, 55, 66};
        int S = 100;
        ArrayList<BigDecimal[]> C = new ArrayList<>();
        ArrayList<BigDecimal[]> b = new ArrayList<>();
        BigDecimal[] c = new BigDecimal[a.length+1];
        Arrays.fill(c,zero);

        //step 1
        double m = ceil(sqrt(a.length)*0.5);

        //step 2
        for (int i = 0; i < a.length; i++)
        {
            c[i] = new BigDecimal(1);
            c[a.length] = new BigDecimal(m*a[i]);
            C.add(c.clone());
            Arrays.fill(c,zero);
        }
        Arrays.fill(c,half);
        c[a.length] = new BigDecimal(S*m);
        C.add(c.clone());

        LLL lll = new LLL(C);
        b = lll.getResult();

        //step 3
    }

}
