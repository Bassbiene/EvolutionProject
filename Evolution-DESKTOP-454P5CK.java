package alg;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class evolution {

	public static final int POPULATIONSGROESSE = 10;
	public static final int WERTEBEREICH_VON = -500;
	public static final int WERTEBEREICH_BIS = 500;
	
	public static double inputFunction(double x1, double x2, double x3, double x4){
		
		double res = 0;
		
		res = Math.pow((x1 + 10 * x2), 2);
		res += 5 * Math.pow((x3 - x4), 2);
		res += Math.pow((x2 - 2* x3), 4);
		res += 10 * Math.pow((x1 - x4),4);
		
		return res;
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		evolve(4);
	}
	
	/***
	 * 
	 * @param n Anzahl der Variablen der inputFunction
	 */
	public static void evolve(int n){
		// TODO: schön machen - Schnittstelle für InputFunktionsklassen
		
		// Schritt 2 - wir Starten mit Generation 0
		int t = 0;
		
		List<double[]> population = new ArrayList();
		// Schritt 3 - Erzeugen der Population
		
		Random r = new Random();
		
		for (int i = 0; i < evolution.POPULATIONSGROESSE; i++){
			
			double[] individuum = new double[n];
			for (int j = 0; j < n; j++){
				//ein Allel eines Individuums mit einem Wert im Wertebereich belegen
				individuum[j] = (double) r.nextInt((evolution.WERTEBEREICH_BIS - evolution.WERTEBEREICH_VON + 1) + evolution.WERTEBEREICH_VON);				
			}
			population.add(individuum);
		}
		System.out.println("test");
		
	}

	public static void printPopulation(List<double[]> population){
		for (double[] individuum : population){
			
		}
	}
	
}
