package alg;

import java.nio.channels.FileLockInterruptionException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class Evolution {

	public static final int POPULATIONSGROESSE = 10;
	public static final int ANZAHL_KINDER = 40;
	public static final int WERTEBEREICH_VON = -500;
	public static final int WERTEBEREICH_BIS = 500;
	
	/***
	 * Die Eingabe-Funktion, für ein Problem, das gelöst werden soll.
	 * @param gen
	 * @return
	 */
	public static double inputFunction(double[] gen){
		
		double res = 0;
		
		res = Math.pow((gen[0] + 10 * gen[1]), 2);
		res += 5 * Math.pow((gen[2] - gen[3]), 2);
		res += Math.pow((gen[1] - 2* gen[2]), 4);
		res += 10 * Math.pow((gen[0] - gen[3]),4);
		
		return res;
	}
	
	/**
	 * Funktion zum Testen einer Lösung
	 * (Bei Minimalwertproblem = Inputfunktion)
	 * @param gen
	 * @return
	 */
	public static double fitnessFunction(double[] gen){
		
		return inputFunction(gen);
	}
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub

		evolve(4);
	}
	
	/***
	 * 
	 * @param n Anzahl der Variablen der inputFunction
	 * @throws Exception 
	 */
	public static void evolve(int n) throws Exception{
		// TODO: schön machen - Schnittstelle für InputFunktionsklassen
		
		// Schritt 2 - wir Starten mit Generation 0
		int t = 0;
		
		List<Individuum> population = new ArrayList();
		// Schritt 3 - Erzeugen der Population
		
		Random r = new Random();
		
		for (int i = 0; i < Evolution.POPULATIONSGROESSE; i++){
			
			double[] allele = new double[n];
			for (int j = 0; j < n; j++){
				//ein Allel eines Individuums mit einem Wert im Wertebereich belegen
				allele[j] = (double) r.nextInt(Evolution.WERTEBEREICH_BIS - Evolution.WERTEBEREICH_VON + 1) + Evolution.WERTEBEREICH_VON;				
			}
			population.add(new Individuum(allele));
		}

		// Ausgabe der Urpopulation
		System.out.println("Urpopulation:");
		printPopulation(population);
		
		// Schritt 4 - Bewertung der Urpopulation
		for (Individuum individuum : population){
			// Bei Nullstellen-Problem müsste man jetzt zuerst eine Lösung aus der Eingabefunktion holen
			// und diese dann mit der Fitnessfunktion testen
			individuum.setFitness(fitnessFunction(individuum.getGenom()));
		}
		
		// Schritt 5 - Leere Kindgeneration anlegen
		List<Individuum> kindgeneration = new ArrayList();
		
		// Sortieren der Eltern nach Fitness
		population.sort(new FitnessComparator());
		System.out.println("");
		printPopulation(population);
			
		// durchgehen der Liste vom besten zum schlechtesten
		// Berechnen der Wahrsch. für jedes Individuum 
		// rangbasierte Selektion
		for (int i = 0; i < Evolution.POPULATIONSGROESSE; i++){
			
			Individuum currentIndividuum = population.get(i);
			currentIndividuum.setWahrscheinlichkeit(berechneWahrscheinlichkeit(i+1, Evolution.POPULATIONSGROESSE));
			if (i > 0){
				
				currentIndividuum.setWahrsch_von(population.get(i-1).getWahrsch_bis());
				currentIndividuum.setWahrsch_bis(currentIndividuum.getWahrsch_von() + currentIndividuum.getWahrscheinlichkeit());				
			}else{
				// für das erste Element
				currentIndividuum.setWahrsch_von(0.0);
				currentIndividuum.setWahrsch_bis(currentIndividuum.getWahrscheinlichkeit());
			}
		}
		System.out.println("");
		printPopulation(population);
		
		// Schritt 6 - 12 Iteration
		for (int i = 1; i<= Evolution.ANZAHL_KINDER; i++){
			
			// Schritt 7 - Eltern auswähelen
			// Elternselektion => Zufällig			
			// Zufallsselektion - rangbasierte Seleketion für Minimierungsprobleme
			
			Individuum neuesKind = null;
			
			while(neuesKind == null){
				
				Individuum elternteil_A = selektiereElternteil(population);
				Individuum elternteil_B = selektiereElternteil(population);
				
				// ToDo: Eltern dürfen nicht beide das gleiche Individuum sein
				
				// Schritt 8 - ein Kind durch Rekombination erzeugen
				
				// Schritt 9 - Kind gegebenenfalls mutieren
				
				// Schritt 10 - neues Kind der Kindergeneration hinzufügen
			}		
			
			// Schritt 13 - Bewerte alle Individuuen der Kindergeneration
			
			// Schritt 14 - Übergang auf die Nachfolgegeneration
			
		}
	}

	public static Individuum selektiereElternteil(List<Individuum> population) throws Exception{
		
		// Erzeugen einer Zufallszahl zwischen 0 und 1.0
		// => muss so kompliziert gemacht werden, weil bei random die 1.0 exclusive ist		
		Random r = new Random();
		double z = (double)r.nextInt(1000000001) / 1000000000;
		
		for (Individuum currentIndividuum : population){
			// Roulette Auswahlverfahren
			if(currentIndividuum.getWahrsch_von() <= z || z < currentIndividuum.getWahrsch_bis()){
				return currentIndividuum;
				// currentIndividuum ist als Elternteil gewählt
			}
		}
		
		throw new Exception("Es wurde kein Individuum mit der gewürfelten Zufallszahl gefunden. Dies kann nicht sein.");
	}
	
	public static void printPopulation(List<Individuum> population){
		for (Individuum individuum : population){
			String currentIndividuum = "";
			
			for (int i = 0; i < individuum.getGenom().length; i++){
				currentIndividuum += individuum.getGenom()[i] + "  ";				
			}
			currentIndividuum += " - Fitness: " + individuum.getFitness();
			currentIndividuum += " - Wahrsch.: " + individuum.getWahrscheinlichkeit();
			currentIndividuum += " - Wahrsch_von: " + individuum.getWahrsch_von();
			currentIndividuum += " - Wahrsch_bis: " + individuum.getWahrsch_bis();
			
			System.out.println(currentIndividuum);
		}
	}
	
	/***
	 * Berechnung der Wahrsch. für rangbasierte Selektion für Minimierungsporbleme
	 * @param i Rang des Individuums, dessen Wahrsch. berechnet werden soll
	 * @param r Anzahl der Individuuen der Population
	 * @return
	 */
	public static Double berechneWahrscheinlichkeit(double i, double r){
		double res = ((2 / r) * (1 - ((i-1)/(r-1))));
		
		return Double.valueOf(res);
	}
	
}
