using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using System.Collections.Concurrent;
using ExcelDna.Integration;

namespace TP_VBA_.NET
{
    static class Constants
    {
        public const int panier = 3;
    }

    public class Option
    {
        #region  Constructor

        public Option()
        {
            option_type = "vanille";
            price_option = 100.0;
            maturity = 1.0;
            strike = 110.0;
            volatility = 30.0 / 100;
            risk_rate = 5.0 / 100;
            call_value = 0.0;
            put_value = 0.0;
            payoff = 0.0;
        }

        public Option(double S, double T, double K, double sigma, double r)
        {
            option_type = "vanille";
            price_option = S;
            maturity = T;
            strike = K;
            volatility = sigma / 100;
            risk_rate = r / 100;
            call_value = 0.0;
            put_value = 0.0;
            payoff = 0.0;
        }

        #endregion

        #region Static Properties

        public static string option_type { get; set; }
        public static double price_option { get; set; }
        public static double maturity { get; set; }
        public static double strike { get; set; }
        public static double volatility { get; set; }
        public static double risk_rate { get; set; }
        public static double call_value { get; set; }
        public static double put_value { get; set; }
        public static double payoff { get; set; }

        #endregion

        #region Static Methods
		//Call of Black Scholes
        public static double CallValue(double S, double K,double T, double sigma, double r)
        {
            price_option = S;
            maturity = T;
            strike = K;
            volatility = sigma / 100;
            risk_rate = r / 100;
            var normal = new MathNet.Numerics.Distributions.Normal(0.0, 1.0);
            var d1 = ((Math.Log(price_option / strike) + maturity * (risk_rate + ((volatility * volatility)/ 2)))) / (volatility * Math.Sqrt(maturity));
            var d2 = d1 - volatility * Math.Sqrt(maturity);
            call_value = price_option * normal.CumulativeDistribution(d1) - strike * Math.Exp(-risk_rate * maturity) * normal.CumulativeDistribution(d2);
            return call_value;
        }
        #endregion

        #region Methods
        public double GeneratePayoffCall(double S, double K)
        {
            return Math.Max(0, S - K);
        }

        #endregion
    }
	
	//Generator of Normal law
    public class BoxMullerGenerator
    {
        #region Fields
            private double? m_already_generated;
            private readonly Random rand;
        #endregion

        #region Constructors
            public BoxMullerGenerator()
            {
                rand = new Random();
            }
			
			public BoxMullerGenerator(int seed)
            {
                rand = new Random(seed);
            }
        #endregion
	
		#region Property

		public double NextDouble
		{
			get { return generate(); }
		}
		#endregion

        #region Public Methods
            public double generate()
            {
			  //If value is not null so we return the value saved
                if (m_already_generated.HasValue)
                {
                    var res = m_already_generated.Value;
                    m_already_generated = null;
                    return res;
                }
			  //generation of tow value
                var u1 = rand.NextDouble();
                var u2 = rand.NextDouble();
			  
                var phase = 2.0 * Math.PI * u1;
                var module = Math.Sqrt(-2 * Math.Log(u2));
              
			  //BoxMuller Formula
                var c = module * Math.Cos(phase);
                var s = module * Math.Sin(phase);
			  //save of second value
                m_already_generated = c;

                return s;
            }
        #endregion
    }
	//Scholesky Brownian Movement
    public class Scholesky
    {
        #region Fields
            private BoxMullerGenerator normal_generator;
            private Matrix<double> matrix_covariance;
        #endregion

        #region Constructors
            public Scholesky()
            {
                normal_generator = new BoxMullerGenerator();
                matrix_covariance = Matrix<double>.Build.Random(Constants.panier, Constants.panier);                               
            }
        #endregion

        #region Public Methods
		// cov(X1,X2)=E[X1X2]-E[X1]E[X2]  -> correlation (X1,X2)= cov(X1,X2)/Racine(Var(X1))*Racine(Var(X2))  
            public Vector<double> generate(double correlation01, double correlation02, double correlation12, double volatility0, double volatility1, double volatility2)
            {
                var rand_vector=Vector<double>.Build.Dense(Constants.panier);
                var scholesky_matrix = Matrix<double>.Build.Random(Constants.panier, Constants.panier);
                var multi_var_gaussian = Vector<double>.Build.Dense(Constants.panier);

                for (int i = 0; i < Constants.panier; ++i)
                {
                    rand_vector.At(i,normal_generator.generate());
                }

                matrix_covariance.At(0, 0, volatility0 * volatility0);
                matrix_covariance.At(1, 1, volatility1 * volatility1);
                matrix_covariance.At(2, 2, volatility2 * volatility2);
                matrix_covariance.At(0, 1, correlation01 * volatility0 * volatility1);
                matrix_covariance.At(0, 2, correlation02 * volatility0 * volatility2);
                matrix_covariance.At(1, 0, correlation01 * volatility0 * volatility1);
                matrix_covariance.At(1, 2, correlation12 * volatility1 * volatility2);
                matrix_covariance.At(2, 0, correlation02 * volatility0 * volatility2);
                matrix_covariance.At(2, 1, correlation12 * volatility1 * volatility2);

                var cholesky = matrix_covariance.Cholesky();
                scholesky_matrix=cholesky.Solve(matrix_covariance);

                multi_var_gaussian = scholesky_matrix.Multiply(rand_vector);
           
                return multi_var_gaussian;             
            }
        #endregion
    }

    public class Pricer
    {
		//Monte Carlo
        public static ConcurrentBag<double>[] generate_monte_carlo(double[] array_spot, double[] exp_drift, Vector<double>[] alea, double[] sqrt_maturity_times_volatility, int nb_iteration, int nb_option)
        {
            var results = new ConcurrentBag<double>[nb_option];    
            for (int j = 0; j < nb_option; ++j)
            {
                results[j] = new ConcurrentBag<double>();
                Parallel.For(0, nb_iteration, i =>
                {
					//price of action at T time (theorem of MonteCarlo)
                    var S = array_spot[j] * Math.Exp(exp_drift[j] + alea[i].At(j) * sqrt_maturity_times_volatility[j]);
                    results[j].Add(S);
                });
            }
            return results;
        }
		//Pricer Call
        public static object generate_pricer_call(double spot0, double strike, double maturity0,double volatility0,  double risk_rate0, int nb_iteration)
        {
            return ExcelAsyncUtil.Run("WebSnippetAsync", new object[] { spot0, strike, maturity0, volatility0, risk_rate0, nb_iteration },
			delegate
			{
				var baskets_of_assets = new Option(spot0, strike, maturity0, volatility0, risk_rate0);
				var normal_generator2 = new BoxMullerGenerator();
				var rand = 0.0;

				var array_spot = new double[] { spot0 };
				var sqrt_maturity_times_volatility = new double[1];
				var exp_drift = new double[1];
				var alea = new Vector<double>[nb_iteration + 1];

				var results = new ConcurrentBag<double>[1];
				var payoff = new List<double>();
				var expectation = new double[1];


				exp_drift[0] = (risk_rate0 / 100.0 - (volatility0 / 100.0 * volatility0 / 100.0) / 2.0) * maturity0;
				sqrt_maturity_times_volatility[0] = Math.Sqrt(maturity0) * (volatility0 / 100.0);

				for (int i = 0; i < nb_iteration + 1; ++i)
				{
					alea[i] = Vector<double>.Build.Dense(1);
					rand = normal_generator2.generate();
					alea[i].At(0, rand);
				}

				results = generate_monte_carlo(array_spot, exp_drift, alea, sqrt_maturity_times_volatility, nb_iteration, 1);

				foreach (double element in results[0])
				{
					payoff.Add(Math.Max(0, element - strike));
				}

				expectation[0] = payoff.Average() * Math.Exp(-risk_rate0 * maturity0 / 100.0);

				return expectation[0];
			});
        }
		//Pricer Put
        public static object generate_pricer_basket(double spot0, double spot1, double spot2, double maturity0, double maturity1, double maturity2, double volatility0, double volatility1, double volatility2, double risk_rate0, double risk_rate1, double risk_rate2, double weight1, double weight2, double weight3, double correlation01, double correlation02, double correlation12, int nb_iteration)
        {
            return ExcelAsyncUtil.Run("WebSnippetAsync", new object[] {  spot0,  spot1,  spot2,  maturity0,  maturity1,  maturity2, volatility0,  volatility1,  volatility2,  risk_rate0,  risk_rate1,  risk_rate2,  weight1, weight2,  weight3, correlation01, correlation02, correlation12, nb_iteration },
			delegate
			{
				var normal_generator = new Scholesky();

				var array_maturity = new double[] { maturity0, maturity1, maturity2 };
				var array_volatility = new double[] { volatility0, volatility1, volatility2 };
				var array_risk_rate = new double[] { risk_rate0, risk_rate1, risk_rate2 };

				var array_spot = new double[] { spot0, spot1, spot2 };
				var sqrt_maturity_times_volatility = new double[Constants.panier];
				var exp_drift = new double[Constants.panier];
				var alea = new Vector<double>[nb_iteration + 1];

				var results = new ConcurrentBag<double>[Constants.panier];
				var results_square = new List<double>[Constants.panier];
				var spot_T = new double[Constants.panier];
				var standard_deviation = new double[Constants.panier];
				var global_spot_T = 0.0;
				var global_standard_deviation = 0.0;
				var global = new double[2];

				for (int i = 0; i < Constants.panier; ++i)
				{
					exp_drift[i] = (array_risk_rate[i] / 100 - (array_volatility[i] / 100 * array_volatility[i] / 100) / 2.0) * array_maturity[i];
					sqrt_maturity_times_volatility[i] = Math.Sqrt(array_maturity[i]) * (array_volatility[i] / 100);
					results[i] = new ConcurrentBag<double>();
					results_square[i] = new List<double>();
				}

				for (int i = 0; i < nb_iteration + 1; ++i)
				{
					alea[i] = normal_generator.generate(correlation01, correlation02, correlation12, volatility0, volatility1, volatility2);
				}

				results = generate_monte_carlo(array_spot, exp_drift, alea, sqrt_maturity_times_volatility, nb_iteration, Constants.panier);

				for (int j = 0; j < Constants.panier; ++j)
				{
					foreach (double element in results[j])
					{
						results_square[j].Add(element * element);
					}

					spot_T[j] = results[j].Average();
					standard_deviation[j] = Math.Sqrt(results_square[j].Average() - (spot_T[j] * spot_T[j]));
				}

				global_spot_T = spot_T[0] * weight1 + spot_T[1] * weight2 + spot_T[2] * weight3;
				global_standard_deviation = standard_deviation[0] * weight1 + standard_deviation[1] * weight2 + standard_deviation[2] * weight3;

				global[0] = global_spot_T;
				global[1] = global_standard_deviation;
				return global;
			});
        }       
    }       
}
