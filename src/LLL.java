import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;

import static java.lang.Math.max;

public class LLL {

    private static ArrayList<BigDecimal[]> b;
    private static ArrayList<BigDecimal[]> b_new = new ArrayList<>();
    private static BigDecimal[] B;
    private static BigDecimal half = new BigDecimal("0.5");
    private static BigDecimal half_h = new BigDecimal("0.75");
    private static BigDecimal zero = new BigDecimal("0");
    private static BigDecimal two = new BigDecimal("2");
    private static MathContext mc = new MathContext(5);

    private static void Gram_Schmidt()
    {
        B = new BigDecimal[b.size()];
        b_new.add(0,b.get(0));
        B[0] = norm(b_new.get(0)).pow(2,mc);
        for(int i=1; i < b.size(); i++)
        {
            b_new.add(i,b.get(i));
            for(int j=0; j < i; j++)
            {
                BigDecimal mu = scalar_multiply(b.get(i), b_new.get(j)).divide(B[j], 5, BigDecimal.ROUND_HALF_UP);
                BigDecimal[] b_temp = subtraction(b_new.get(i), scalar_multiply(mu, b_new.get(j)));
                b_new.set(i,b_temp);
            }
            B[i] = norm(b_new.get(i)).pow(2,mc);
        }
    }

    private static BigDecimal norm(BigDecimal[] vector)
    {
        return sqrt(scalar_multiply(vector, vector));
    }

    private static BigDecimal[] scalar_multiply(BigDecimal c, BigDecimal[] b)
    {
        BigDecimal[] result = new BigDecimal[b.length];
        for(int i=0; i < b.length; i++)
            result[i] = c.multiply(b[i],mc);
        return result;
    }

    private static BigDecimal scalar_multiply(BigDecimal[] a, BigDecimal[] b)
    {
        BigDecimal result = new BigDecimal("0");
        for(int i=0; i < a.length; i++)
            result = result.add(a[i].multiply(b[i],mc));
        return result;
    }

    private static BigDecimal[] subtraction(BigDecimal[] a, BigDecimal[] b)
    {
        BigDecimal[] result = new BigDecimal[a.length];
        for(int i=0; i < a.length; i++)
            result[i] = a[i].subtract(b[i]);
        return result;
    }

    private static BigDecimal sqrt(BigDecimal value) {
        BigDecimal x = new BigDecimal(Math.sqrt(value.doubleValue()));
        return x.add(new BigDecimal(value.subtract(x.multiply(x,mc)).doubleValue() / (x.doubleValue() * 2.0)));
    }

    private static void LLL_algorithm()
    {
        //step 1
        Gram_Schmidt();

        //debug
        //step 2
        int k = 1;
        BigDecimal r;
        BigDecimal mu_temp;
        BigDecimal[][] mu = new BigDecimal[b.size()][b.size()];

        for (int i = 0; i < b.size(); i++)
            for (int j = 0; j < b.size(); j++)
                mu[i][j] = scalar_multiply(b.get(i), b_new.get(j)).divide(B[j], 5, BigDecimal.ROUND_HALF_UP);

        while(true)
        {
            //step 3
            if (mu[k][k - 1].abs().compareTo(half) == 1)
            {
                //step 3.1
                r = (mu[k][k - 1].compareTo(zero) == 1) ? mu[k][k - 1].add(half).setScale(0, RoundingMode.FLOOR) : half.subtract(mu[k][k - 1]).setScale(0, RoundingMode.FLOOR).negate();
                //step 3.2
                b.set(k, subtraction(b.get(k), scalar_multiply(r, b.get(k - 1))));
                //step 3.3
                for (int j = 0; j < k - 1; j++)
                    mu[k][j] = mu[k][j].subtract(r.multiply(mu[k - 1][j],mc));
                //step 3.4
                mu[k][k - 1] = mu[k][k - 1].subtract(r);
            }

            //step 4
            if (B[k].compareTo(half_h.subtract(mu[k][k - 1].pow(2,mc)).multiply(B[k - 1],mc)) == -1) {
                //step 4.1
                mu_temp = mu[k][k - 1];
                BigDecimal B_temp = B[k].add(mu_temp.pow(2,mc).multiply(B[k - 1],mc));
                mu[k][k - 1] = mu_temp.multiply(B[k - 1],mc).divide(B_temp, 5, BigDecimal.ROUND_HALF_UP);
                B[k] = B[k - 1].multiply(B[k],mc).divide(B_temp, 5, BigDecimal.ROUND_HALF_UP);
                B[k - 1] = B_temp;
                //step 4.2
                BigDecimal[] temp = b.get(k);
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
                    BigDecimal t = mu[s][k];
                    mu[s][k] = mu[s][k - 1].subtract(mu_temp.multiply(t,mc));
                    mu[s][k - 1] = t.add(mu[k][k - 1].multiply(mu[s][k],mc));
                }
                //step 4.5
                k = max(1, k - 1);
            }
            //step 5
            else
            {
                //step 5.1
                for (int l = k-1; l >= 0; l--)
                    if(mu[k][l].abs().compareTo(half) == 1)
                    {
                        //step 5.1.1
                        r = (mu[k][l].compareTo(zero) == 1) ? mu[k][l].add(half).setScale(0, RoundingMode.FLOOR) : half.subtract(mu[k][l]).setScale(0, RoundingMode.FLOOR).negate();
                        //step 5.1.2
                        b.set(k, subtraction(b.get(k), scalar_multiply(r, b.get(l))));
                        //step 5.1.3
                        for(int j = 0; j < l; j++)
                            mu[k][j] = mu[k][j].subtract(r.multiply(mu[l][j],mc));
                        //step 5.1.4
                        mu[k][l] = mu[k][l].subtract(r);
                    }
                //step 5.2
                k++;
            }

            if(k > b.size() - 1)
                break;
        }

        //LLL result
        BigDecimal modulo = new BigDecimal(b.size());
        for (BigDecimal[] aB : b) {
            for (BigDecimal anAB : aB) System.out.print(anAB + " ");
            System.out.println();
        }
    }

    ArrayList<BigDecimal[]> getResult()
    {
        return b;
    }

    public LLL(ArrayList<BigDecimal[]> vectors)
    {
/*        BigDecimal[] b1 = {BigDecimal.ONE,BigDecimal.ZERO,BigDecimal.ONE,two};
        BigDecimal[] b2 = {BigDecimal.ONE,BigDecimal.ONE.negate(),two,zero};
        BigDecimal[] b3 = {BigDecimal.ONE.negate(),two,zero,BigDecimal.ONE};
        b = new ArrayList<>();
        b.add(b1);
        b.add(b2);
        b.add(b3);*/
        b = new ArrayList<>(vectors);
        LLL_algorithm();
    }
}
