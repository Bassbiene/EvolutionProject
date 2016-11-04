package alg;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Evolution {

	public static final int MAX_GENERATIONEN = 2000;
	public static final int POPULATIONSGROESSE = 100;
	public static final int ANZAHL_KINDER = 100;
	public static final int WERTEBEREICH_VON = -500;
	public static final int WERTEBEREICH_BIS = 500;
	public static final double REKOMBINATIONSWAHRSCH = 0.3;
	public static final double MUTATIONSWAHRSCH = 0.1;
	public static final Elternselektion ELTERNSELEKTION = Elternselektion.ZufaelligGleichverteilt;
	public static final Rekombination REKOMBINATION = Rekombination.Intermediaer;
	
	/***
	 * Die Eingabe-Funktion, fuer ein Problem, das geloest werden soll.
	 * 
	 * @param gen
	 * @return
	 */
	public static double inputFunction(double[] gen) {

		double res = 0;

		res = Math.pow((gen[0] + 10 * gen[1]), 2);
		res += 5 * Math.pow((gen[2] - gen[3]), 2);
		res += Math.pow((gen[1] - 2 * gen[2]), 4);
		res += 10 * Math.pow((gen[0] - gen[3]), 4);

		return res;
	}

	/**
	 * Funktion zum Testen einer Lï¿½sung (Bei Minimalwertproblem =
	 * Inputfunktion)
	 * 
	 * @param gen
	 * @return
	 */
	public static double fitnessFunction(double[] gen) {

		return inputFunction(gen);
	}

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub

		evolve(4);
	}

	/***
	 * 
	 * @param n
	 *            Anzahl der Variablen der inputFunction (= Groesse des Genoms)
	 * @throws Exception
	 */
	public static void evolve(int n) throws Exception {
		// TODO: schï¿½n machen - Schnittstelle fï¿½r InputFunktionsklassen

		List<Individuum> population = new ArrayList<Individuum>();

		// Schritt 3 - Erzeugen der Urpopulation
		Random r = new Random();

		for (int i = 0; i < Evolution.POPULATIONSGROESSE; i++) {

			double[] allele = new double[n];
			for (int j = 0; j < n; j++) {
				// ein Allel eines Individuums mit einem Wert im
				// Wertebereich belegen
				allele[j] = (double) r.nextInt(Evolution.WERTEBEREICH_BIS - Evolution.WERTEBEREICH_VON + 1)
						+ Evolution.WERTEBEREICH_VON;
			}
			population.add(new Individuum(allele));
		}

		// Ausgabe der Urpopulation
		// System.out.println("Urpopulation:");
		// printPopulation(population);

		// Schritt 2 - wir Starten mit Generation 0
		for (int t = 0; t < Evolution.MAX_GENERATIONEN; t++) {

			System.out.println("Generation: " + t);

			// Schritt 4 - Bewertung der Population
			printCurrentStep("Schritt 4 - Bewertung der Population");
			bewertePopulation(population);

			// Schritt 5 - Leere Kindgeneration anlegen
			printCurrentStep("Schritt 5 - Leere Kindgeneration anlegen");
			List<Individuum> kindgeneration = new ArrayList<>();

			// Sortiert die Individuuen nach Fitness
			// und ordnet Ihnen die rangbasierte Wahrscheinlichkeit für die
			// Rouletteselektion zu
			if (Evolution.ELTERNSELEKTION == Elternselektion.Rouletteverfahren) {
				bereiteRouletteSelektionVor(population);
			}

			// Schritt 6 - 12 Iteration
			for (int i = 1; i <= Evolution.ANZAHL_KINDER; i++) {

				// Schritt 7 - Eltern auswï¿½helen
				// Elternselektion => Zufï¿½llig
				// Zufallsselektion - rangbasierte Seleketion fï¿½r
				// Minimierungsprobleme

				printCurrentStep("Schritt 7 - Eltern auswaehlen");
				Individuum neuesKind = null;

				while (neuesKind == null) {

					Individuum elternteil_A = null;
					Individuum elternteil_B = null;

					switch (Evolution.ELTERNSELEKTION) {
					case Rouletteverfahren:
						elternteil_A = selektiereIndividuumRoulette(population);
						elternteil_B = selektiereIndividuumRoulette(population);
						// Eltern duerfen nicht beide das gleiche Individuum sein
						while (elternteil_A == elternteil_B) {
							elternteil_B = selektiereIndividuumRoulette(population);
						}
						break;
					case ZufaelligGleichverteilt:
						elternteil_A = selektiereIndividuumZufaelligGleichverteilt(population);
						elternteil_B = selektiereIndividuumZufaelligGleichverteilt(population);
						// Eltern duerfen nicht beide das gleiche Individuum sein
						while (elternteil_A == elternteil_B) {
							elternteil_B = selektiereIndividuumZufaelligGleichverteilt(population);
						}
						break;
					}

					// die Wahrsch. sorgt dafür, dass sich nicht nur die besten
					// x rekombinieren, sond. auch andere eine Chance haben
					double zufallRekombination = r.nextDouble();
					if (zufallRekombination < Evolution.REKOMBINATIONSWAHRSCH) {
						// Schritt 8 - ein Kind durch Rekombination erzeugen
						printCurrentStep("Schritt 8 - ein Kind durch Rekombination erzeugen");
						
						switch (Evolution.REKOMBINATION){
						case Arithmetisch:
							neuesKind = rekombiniereElternArithmetisch(elternteil_A, elternteil_B);
							break;
						case Intermediaer:
							neuesKind = rekombiniereElternIntermediaer(elternteil_A, elternteil_B);
							break;
						}
						
					}

				}

				// Schritt 9 - Kind gegebenenfalls mutieren
				printCurrentStep("Schritt 9 - Kind gegebenenfalls mutieren");
				double zufallMutation = r.nextDouble();
				if (zufallMutation < Evolution.MUTATIONSWAHRSCH) {

					mutiereIndividuum(neuesKind);
				}

				// Schritt 10 - neues Kind der Kindergeneration hinzufuegen
				printCurrentStep("Schritt 10 - neues Kind der Kindergeneration hinzufuegen");
				kindgeneration.add(neuesKind);

			}
			// Schritt 13 - Bewerte alle Individuuen der Kindergeneration
			printCurrentStep("Schritt 13 - Bewerte alle Individuuen der Kindergeneration");
			bewertePopulation(kindgeneration);

			// Schritt 14 - Uebergang auf die Nachfolgegeneration
			printCurrentStep("Schritt 14 - Uebergang auf die Nachfolgegeneration");
			List<Individuum> gesamtpopulation = new ArrayList<>();
			// die Gesamtpopulation soll alle Eltern und Kinder enthalten,
			// sodass bei der Umweltselektion auch die Eltern berücksichtigt
			// werden und in die
			// nächste Generation übernommen werden können
			gesamtpopulation.addAll(population);
			gesamtpopulation.addAll(kindgeneration);

			population.clear();

			// Schritt 15 - Selektiere Individuuen für die nächste Generation
			printCurrentStep("Schritt 15 - Selektiere Individuuen für die nächste Generation");
			// Umweltselektion
			bereiteRouletteSelektionVor(gesamtpopulation);
			
			while (population.size() < Evolution.POPULATIONSGROESSE) {

				Individuum selektiertesIndividuum = selektiereIndividuumRoulette(gesamtpopulation);
				if (!population.contains(selektiertesIndividuum)) {
					population.add(selektiertesIndividuum);
				}
			}
		}

		// Schritt 17: Ermittle bestes Individuum nach Ablauf des
		// Evolutionsalgorithmus
		printCurrentStep("Schritt 17: Ermittle bestes Individuum nach Ablauf des Evolutionsalgorithmus");
		bewertePopulation(population);
		population.sort(new FitnessComparator());

		System.out.println("\nDies ist die letzte Population:");
		printPopulation(population);

		System.out.println("\nDies ist das beste Individuum:");
		printIndividuum(population.get(0));

	}

	/**
	 * Sortiert die Individuuen nach Fitness Ordnet Ihnen die rangbasierte
	 * Wahrscheinlichkeit für die Rouletteselektion zu
	 * 
	 * @param population
	 */
	public static void bereiteRouletteSelektionVor(List<Individuum> population) {
		// Sortieren der Eltern nach Fitness
		population.sort(new FitnessComparator());

		// durchgehen der Liste vom besten zum schlechtesten
		// Berechnen der Wahrsch. fï¿½r jedes Individuum
		// rangbasierte Selektion
		for (int i = 0; i < population.size(); i++) {

			Individuum currentIndividuum = population.get(i);
			currentIndividuum.setWahrscheinlichkeit(berechneWahrscheinlichkeit(i + 1, population.size()));
			if (i > 0) {

				currentIndividuum.setWahrsch_von(population.get(i - 1).getWahrsch_bis());
				currentIndividuum
						.setWahrsch_bis(currentIndividuum.getWahrsch_von() + currentIndividuum.getWahrscheinlichkeit());
			} else {
				// fï¿½r das erste Element
				currentIndividuum.setWahrsch_von(0.0);
				currentIndividuum.setWahrsch_bis(currentIndividuum.getWahrscheinlichkeit());
			}
		}
	}

	public static void bewertePopulation(List<Individuum> population) {
		for (Individuum individuum : population) {
			// Bei Nullstellen-Problem mï¿½sste man jetzt zuerst eine Lï¿½sung
			// aus der Eingabefunktion holen
			// und diese dann mit der Fitnessfunktion testen
			individuum.setFitness(fitnessFunction(individuum.getGenom()));
		}
	}

	/**
	 * Waehlt zufaellig ein Gen des Individuums und mutiert dieses innerhalb des
	 * definierten Wertebereichs.
	 * 
	 * @param ind
	 */
	public static void mutiereIndividuum(Individuum ind) {

		double[] genom = ind.getGenom();

		Random r = new Random();

		// wähle zufällig, welches Gen mutiert werden soll
		// wirft eine Zufallszahl zwischen 0 inkl. und genom.length exklusive
		int zufallGen = r.nextInt(genom.length);

		// Erzeuge eine Zufallszahl im definierten Wertebereich
		double zufallsZahl = (double) r.nextInt(Evolution.WERTEBEREICH_BIS - Evolution.WERTEBEREICH_VON + 1)
				+ Evolution.WERTEBEREICH_VON;

		// mutiere das Gen
		genom[zufallGen] = genom[zufallGen] + zufallsZahl;
		// falls der Wertebereich überschritten wird, sorge dafür dass die
		// Schranken eingehalten werden
		if (genom[zufallGen] > Evolution.WERTEBEREICH_BIS) {
			genom[zufallGen] = Evolution.WERTEBEREICH_BIS;
		}
		if (genom[zufallGen] < Evolution.WERTEBEREICH_VON) {
			genom[zufallGen] = Evolution.WERTEBEREICH_VON;
		}

		ind.setGenom(genom);
	}

	/***
	 * Intermediaere Rekombination Rekombiniert Elternteil A und ELternteil B
	 * anhand des arithmetischen mittels ihrer Allele
	 * 
	 * @param elternteil_A
	 * @param elternteil_B
	 * @return
	 */
	public static Individuum rekombiniereElternIntermediaer(Individuum elternteil_A, Individuum elternteil_B) {

		double[] genom_A = elternteil_A.getGenom();
		double[] genom_B = elternteil_B.getGenom();

		double[] genom_Kind = new double[genom_A.length];

		// rekombiniere alle einzelnen Allele (arithmetisches Mittel)
		for (int i = 0; i < genom_A.length; i++) {
			genom_Kind[i] = (genom_A[i] + genom_B[i]) / 2;
		}

		return new Individuum(genom_Kind);
	}

	/**
	 * Arithmetische Rekombination Rekombiniert Elternteil A und ELternteil B
	 * anhand eines zufaelligen Werts zwischen den beieden Allelen von A und B
	 * @param elternteil_A
	 * @param elternteil_B
	 * @return
	 */
	public static Individuum rekombiniereElternArithmetisch(Individuum elternteil_A, Individuum elternteil_B){
		
		double[] genom_A = elternteil_A.getGenom();
		double[] genom_B = elternteil_B.getGenom();

		double[] genom_Kind = new double[genom_A.length];
		
		// rekombiniere alle einzelnen Allele (zufälliger Wert zwischen den beiden Werten)
		Random r = new Random();
		
		for (int i = 0; i < genom_A.length; i++) {
			
			double wertBis = Math.max(genom_A[i] , genom_B[i]);
			double wertVon = Math.min(genom_A[i] , genom_B[i]);
			
			genom_Kind[i] = (double) (r.nextInt((int)(wertBis - wertVon + 1) * 1000)) / 1000 + wertVon;
			
			//genom_Kind[i] = ((double) r.nextInt((int)( genom_B[i] - genom_A[i] + 1) * 1000) + genom_A[i]) / 1000;	
			//allele[j] = (double) r.nextInt(Evolution.WERTEBEREICH_BIS - Evolution.WERTEBEREICH_VON + 1) + Evolution.WERTEBEREICH_VON;
		}
		
		return new Individuum(genom_Kind);
	}
	
	/**
	 * Rouletteverfahren Wählt zufällig ein Individuum aus einer Population
	 * anhand des Rouletteverfahrens
	 * 
	 * @param population
	 *            Die Population, aus denen ein Individuum gewählt werden soll
	 * @return gewähltes Individuum
	 * @throws Exception
	 */
	public static Individuum selektiereIndividuumRoulette(List<Individuum> population) throws Exception {

		// Erzeugen einer Zufallszahl zwischen 0 und 1.0
		// => muss so kompliziert gemacht werden, weil bei random die 1.0
		// exclusive ist
		Random r = new Random();
		double z = (double) r.nextInt(1000000001) / 1000000000;

		for (Individuum currentIndividuum : population) {
			// Roulette Auswahlverfahren
			if (currentIndividuum.getWahrsch_von() <= z && z < currentIndividuum.getWahrsch_bis()) {
				return currentIndividuum;
				// currentIndividuum ist als Elternteil gewï¿½hlt
			}
		}

		throw new Exception(
				"Es wurde kein Individuum mit der gewï¿½rfelten Zufallszahl gefunden. Dies kann nicht sein.");
	}

	/**
	 * Gibt rein zufällig ein Individuum aus der Population zurück.
	 * @param population
	 * @return
	 */
	public static Individuum selektiereIndividuumZufaelligGleichverteilt(List<Individuum> population) {
		
		Random r = new Random();
		int zufallsZahl = r.nextInt(population.size());
		
		return population.get(zufallsZahl);
	}

	/**
	 * Gibt alle Individuuen einer Population auf der Konsole aus.
	 * 
	 * @param population
	 */
	public static void printPopulation(List<Individuum> population) {
		for (Individuum individuum : population) {
			printIndividuum(individuum);
		}
	}

	/**
	 * Gibt das Individuum auf der Konsole aus.
	 * 
	 * @param individuum
	 */
	public static void printIndividuum(Individuum individuum) {
		String ausgabeString = "";

		for (int i = 0; i < individuum.getGenom().length; i++) {
			ausgabeString += individuum.getGenom()[i] + "  ";
		}
		ausgabeString += " - Fitness: " + individuum.getFitness();
		ausgabeString += " - Wahrsch.: " + individuum.getWahrscheinlichkeit();
		ausgabeString += " - Wahrsch_von: " + individuum.getWahrsch_von();
		ausgabeString += " - Wahrsch_bis: " + individuum.getWahrsch_bis();

		System.out.println(ausgabeString);
	}

	/***
	 * Berechnung der Wahrsch. fï¿½r rangbasierte Selektion fï¿½r
	 * Minimierungsporbleme
	 * 
	 * @param i
	 *            Rang des Individuums, dessen Wahrsch. berechnet werden soll
	 * @param r
	 *            Anzahl der Individuuen der Population
	 * @return
	 */
	public static Double berechneWahrscheinlichkeit(double i, double r) {
		double res = ((2 / r) * (1 - ((i - 1) / (r - 1))));

		return Double.valueOf(res);
	}
	
	/**
	 * Gibt den Text des aktuellen Steps aus.
	 * @param text
	 */
	public static void printCurrentStep(String text){
		//System.out.println(text);
	}

}
