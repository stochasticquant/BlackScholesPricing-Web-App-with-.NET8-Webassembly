namespace BlackScholesPricing.Services
{
    public class BlackScholesService
    {
        
        public double CND(double z)
        {
            double p = 0.3275911;
            double a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
            double a4 = -1.453152027, a5 = 1.061405429;

            int sign = z < 0 ? -1 : 1;
            double x = Math.Abs(z) / Math.Sqrt(2.0);
            double t = 1.0 / (1.0 + p * x);
            double erf = 1.0 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x));

            return 0.5 * (1.0 + sign * erf);
        }

        public double BlackScholes(double S, double K, double r, double vol, double T, char optionType)
        {
            double d1 = (Math.Log(S / K) + (r + 0.5 * vol * vol) * T) / (vol * Math.Sqrt(T));
            double d2 = d1 - vol * Math.Sqrt(T);

            if (optionType == 'c')
                return S * CND(d1) - K * Math.Exp(-r * T) * CND(d2);
            else if (optionType == 'p')
                return K * Math.Exp(-r * T) * CND(-d2) - S * CND(-d1);
            else
                throw new ArgumentException("Invalid option type");
        }
    }
}

