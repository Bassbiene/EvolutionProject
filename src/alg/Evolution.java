package alg;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Random;

public class Evolution {

    public static final int GENOMGROESSE = 20;
    public static final int MAX_GENERATIONEN = 2000;
    public static final int POPULATIONSGROESSE = 200;
    public static final int ANZAHL_KINDER = 800;
    public static final int WERTEBEREICH_VON = -512;
    public static final int WERTEBEREICH_BIS = 511;
    public static final double REKOMBINATIONSWAHRSCH = 0.75;
    public static final double MUTATIONSWAHRSCH = 0.1;
    public static final Elternselektion ELTERNSELEKTION = Elternselektion.Rouletteverfahren;
    public static final Umweltselektion UMWELTSELEKTION = Umweltselektion.Deterministisch;
    public static final Rekombination REKOMBINATION = Rekombination.Arithmetisch;
    public static final EingabeFunktion EINGABEFUNKTION = EingabeFunktion.GRIEWANK;
    public static final boolean ZWISCHENERGEBNISSE_AUSGEBEN = false;
    public static final int ABBRUCH_DURCHLAEUFE_OHNE_AENDERUNG = 20;

    public static void main(String[] args) throws Exception {

        // List<Individuum> population = new ArrayList<Individuum>();
        // population.add(new Individuum(new double[] {1.0, 1.0, 1.0, 1.0}));
        // bewertePopulation(population);
        // printPopulation(population);

        evolve(Evolution.GENOMGROESSE);
    }

    /***
     * Die Eingabe-Funktion, fuer ein Problem, das geloest werden soll.
     * 
     * @param genom
     * @return den Funktionswert der Funktion für das gegebene Genom
     */
    public static double inputFunction(double[] genom) {

        double res = 0;

        switch (Evolution.EINGABEFUNKTION) {
        case SKRIPT:
            // Testfunktion aus dem Skript
            res = Math.pow((genom[0] + (10 * genom[1])), 2);
            res += 5 * Math.pow((genom[2] - genom[3]), 2);
            res += Math.pow((genom[1] - (2 * genom[2])), 4);
            res += 10 * Math.pow((genom[0] - genom[3]), 4);
            break;
        case PARABEL:
            // einfache Testfunktion
            res = Math.pow(genom[0] - 2, 2);
            res += 2;
            break;
        case GRIEWANK:
            // Griewank-Funktion
            res = 1 + summiereGriewank(genom);
            res -= multipliziereGriewank(genom);

            break;
        }

        return res;
    }

    /**
     * Hilfsmethode zur Berechnung des zweiten Teils (Produkt) der
     * Griewank-Funktion
     * 
     * @param genom
     * @return
     */
    private static double multipliziereGriewank(double[] genom) {
        double res;

        res = Math.cos((genom[0] / Math.sqrt(1)));

        for (int i = 1; i < genom.length; i++) {
            res *= Math.cos((genom[i] / Math.sqrt(i + 1)));
        }

        return res;
    }

    /**
     * Hilfsmethode zur Berechnung des ersten Teils (Summe) der
     * Griewank-Funktion
     * 
     * @param gen
     * @return
     */
    private static double summiereGriewank(double[] gen) {
        double res = 0;

        for (int i = 0; i < gen.length; i++) {
            double summe = Math.pow(gen[i], 2);
            summe /= (400 * gen.length);
            res += summe;
        }

        return res;
    }

    /**
     * Funktion zum Testen einer L�sung (Bei Minimalwertproblem = Inputfunktion)
     * 
     * @param gen
     * @return
     */
    public static double fitnessFunction(double[] gen) {

        return inputFunction(gen);
    }

    /**
     * 
     * @param n
     *            Anzahl der Variablen der inputFunction (= Groesse des Genoms)
     * @throws Exception
     */
    public static void evolve(int n) throws Exception {
        // TODO: schoen machen - Schnittstelle fuer InputFunktionsklassen

        double championFitness = Double.MAX_VALUE;
        int anzahlDurchlaeufeMitGleichemChampion = -1;

        Random r = new Random();

        // Schritt 3 - Erzeugen der Urpopulation
        List<Individuum> population = new ArrayList<Individuum>();

        erzeugeUrpopulation(n, population);

        // Schritt 4 - Bewertung der Urpopulation
        printCurrentStep("Schritt 4 - Bewertung der Population");
        bewertePopulation(population);

        // Schritt 2 - wir Starten mit Generation 0
        for (int t = 0; t < Evolution.MAX_GENERATIONEN; t++) {

            System.out.println("\n###############################################");
            System.out.println("Generation: " + t);

            // Schritt 5 - Leere Kindgeneration anlegen
            printCurrentStep("Schritt 5 - Leere Kindgeneration anlegen");
            List<Individuum> kindgeneration = new ArrayList<>();

            if (Evolution.ELTERNSELEKTION == Elternselektion.Rouletteverfahren) {
                bereiteRouletteSelektionVor(population);
            }

            if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                // Ausgabe der Urpopulation
                System.out.println("Population " + t + ":");
                printPopulationMitWahrscheinlichkeit(population);
            }

            // Schritt 6 - 12 Iteration
            for (int i = 1; i <= Evolution.ANZAHL_KINDER; i++) {

                // Schritt 7 - Eltern ausw�helen
                // Elternselektion => Zuf�llig
                // Zufallsselektion - rangbasierte Seleketion f�r
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
                        // Eltern duerfen nicht beide das gleiche Individuum
                        // sein
                        while (elternteil_A == elternteil_B) {
                            elternteil_B = selektiereIndividuumRoulette(
                                    population);
                        }
                        break;
                    case ZufaelligGleichverteilt:
                        elternteil_A = selektiereIndividuumZufaelligGleichverteilt(
                                population);
                        elternteil_B = selektiereIndividuumZufaelligGleichverteilt(
                                population);
                        // Eltern duerfen nicht beide das gleiche Individuum
                        // sein
                        while (elternteil_A == elternteil_B) {
                            elternteil_B = selektiereIndividuumZufaelligGleichverteilt(
                                    population);
                        }
                        break;
                    }

                    if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                        // Ausgabe der Urpopulation
                        System.out.println("\n----------------------------");
                        System.out.println("Zwei Eltern wurden gewaehlt: ");
                        printGenomMitWahrscheinlichkeit(elternteil_A);
                        printGenomMitWahrscheinlichkeit(elternteil_B);
                    }

                    // die Wahrsch. sorgt daf�r, dass sich nicht nur die besten
                    // x rekombinieren, sond. auch andere eine Chance haben
                    double zufallRekombination = r.nextDouble();
                    if (zufallRekombination < Evolution.REKOMBINATIONSWAHRSCH) {

                        if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                            System.out.println("\nZufallswert Rekombination: "
                                    + zufallRekombination);
                            System.out.println(
                                    "Die beiden Individuuen werden nun rekombiniert.");
                        }

                        // Schritt 8 - ein Kind durch Rekombination erzeugen
                        printCurrentStep(
                                "Schritt 8 - ein Kind durch Rekombination erzeugen");

                        switch (Evolution.REKOMBINATION) {
                        case Arithmetisch:
                            neuesKind = rekombiniereElternArithmetisch(
                                    elternteil_A, elternteil_B);
                            break;
                        case Intermediaer:
                            neuesKind = rekombiniereElternIntermediaer(
                                    elternteil_A, elternteil_B);
                            break;
                        }

                        if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                            System.out.println(
                                    "\nDie Rekombination ergab folgendes Kind:");
                            printGenom(neuesKind);
                            
                        }

                    } else {
                        if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                            System.out.println("\nZufallswert Rekombination: "
                                    + zufallRekombination);
                            System.out.println(
                                    "Die beiden Individuuen wurden nicht rekombiniert.");
                        }
                    }
                   

                }

                // Schritt 9 - Kind gegebenenfalls mutieren
                printCurrentStep("Schritt 9 - Kind gegebenenfalls mutieren");
                double zufallMutation = r.nextDouble();

                if (zufallMutation < Evolution.MUTATIONSWAHRSCH) {

                    mutiereIndividuum(neuesKind);

                    if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                        System.out.println(
                                "\nZufallszahl Mutation: " + zufallMutation);
                        System.out.println("Das Individuum wurde mutiert:");
                        printGenom(neuesKind);
                    }

                } else {
                    if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                        System.out.println(
                                "\nZufallszahl Mutation: " + zufallMutation);
                        System.out
                                .println("Das Individuum wurde nicht mutiert.");
                    }
                }

                // Schritt 10 - neues Kind der Kindergeneration hinzufuegen
                printCurrentStep(
                        "Schritt 10 - neues Kind der Kindergeneration hinzufuegen");
                kindgeneration.add(neuesKind);

            }
            // Schritt 13 - Bewerte alle Individuuen der Kindergeneration
            printCurrentStep(
                    "Schritt 13 - Bewerte alle Individuuen der Kindergeneration");
            bewertePopulation(kindgeneration);

            // Schritt 14 - Uebergang auf die Nachfolgegeneration
            printCurrentStep(
                    "Schritt 14 - Uebergang auf die Nachfolgegeneration");
            List<Individuum> gesamtpopulation = new ArrayList<>();
            // die Gesamtpopulation soll alle Eltern und Kinder enthalten,
            // sodass bei der Umweltselektion auch die Eltern ber�cksichtigt
            // werden und in die
            // n�chste Generation �bernommen werden k�nnen
            gesamtpopulation.addAll(population);
            gesamtpopulation.addAll(kindgeneration);

            population.clear();

            // Schritt 15 - Selektiere Individuuen f�r die n�chste Generation
            printCurrentStep(
                    "Schritt 15 - Selektiere Individuuen f�r die n�chste Generation");

            // Umweltselektion
            switch (Evolution.UMWELTSELEKTION) {
            case Roulette:

                bereiteRouletteSelektionVor(gesamtpopulation);
                if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                    System.out.println("\nUmweltselektion: Rouletteverfahren");
                    System.out.println(
                            "Aktuelle Gesamtpopulation mit Selektionswahrscheinlichkeiten:");
                    printPopulation(gesamtpopulation);
                }

                while (population.size() < Evolution.POPULATIONSGROESSE) {

                    Individuum selektiertesIndividuum = selektiereIndividuumRoulette(
                            gesamtpopulation);
                    if (!population.contains(selektiertesIndividuum)) {
                        population.add(selektiertesIndividuum);
                    }
                }
                break;

            case Deterministisch:
                // Sortiere die Gesamtpopulation entsprechend ihrer Fitness und
                // w�hle nur die besten
                // fuer die naechste Generation
                gesamtpopulation.sort(new FitnessComparator());

                if (Evolution.ZWISCHENERGEBNISSE_AUSGEBEN) {
                    System.out.println("\nUmweltselektion: Deterministisch");
                    System.out
                            .println("Aktuelle Gesamtpopulation mit Fitness:");
                    printPopulation(gesamtpopulation);
                }

                for (int i = 0; i < Evolution.POPULATIONSGROESSE; i++) {
                    population.add(gesamtpopulation.get(i));
                }

                break;
            }

            // Pr�fe, wie oft sich die Fitness des Champions ver�ndert hat
            // Beende Evolutionsalgorithmus, wenn sich die Fitness entspr. der
            // Abbruchanzahl der Versuche
            // nicht mehr veraendert hat
            population.sort(new FitnessComparator());
            if (anzahlDurchlaeufeMitGleichemChampion == -1) {
                championFitness = population.get(0).getFitness();
                anzahlDurchlaeufeMitGleichemChampion = 0;

            } else if (anzahlDurchlaeufeMitGleichemChampion < Evolution.ABBRUCH_DURCHLAEUFE_OHNE_AENDERUNG) {

                double newChampionFitness = population.get(0).getFitness();
                if (championFitness == newChampionFitness) {
                    anzahlDurchlaeufeMitGleichemChampion += 1;
                } else {
                    championFitness = newChampionFitness;
                    anzahlDurchlaeufeMitGleichemChampion = 0;
                }
            } else {
                System.out.println(
                        "\nAnzahl durchlaeufe fuer den besten Champion: " + (t
                                - Evolution.ABBRUCH_DURCHLAEUFE_OHNE_AENDERUNG));
                break;
            }

        }

        // Schritt 17: Ermittle bestes Individuum nach Ablauf des
        // Evolutionsalgorithmus
        printCurrentStep(
                "Schritt 17: Ermittle bestes Individuum nach Ablauf des Evolutionsalgorithmus");

        // System.out.println("\nDies ist die letzte Population:");
        // printPopulation(population);

        System.out.println("\nDies ist das beste Individuum:");
        printIndividuum(population.get(0));

    }

    

    /**
     * Erzeugt eine Population mit n Individuuen
     * 
     * @param n
     * @param population
     */
    public static void erzeugeUrpopulation(int n, List<Individuum> population) {

        Random r = new Random();

        for (int i = 0; i < Evolution.POPULATIONSGROESSE; i++) {

            double[] allele = new double[n];
            for (int j = 0; j < n; j++) {
                // ein Allel eines Individuums mit einem Wert im
                // Wertebereich belegen
                double zufallszahl = r.nextDouble();
                allele[j] = zufallszahl * (Evolution.WERTEBEREICH_BIS
                        - Evolution.WERTEBEREICH_VON) + Evolution.WERTEBEREICH_VON;
            }
            
            
            population.add(new Individuum(allele));
        }
    }

    /**
     * Sortiert die Individuuen nach Fitness Ordnet Ihnen die rangbasierte
     * Wahrscheinlichkeit f�r die Rouletteselektion zu
     * 
     * @param population
     */
    public static void bereiteRouletteSelektionVor(
            List<Individuum> population) {
        // Sortieren der Individuen nach Fitness
        population.sort(new FitnessComparator());

        // durchgehen der Liste vom besten zum schlechtesten
        // Berechnen der Wahrsch. fuer jedes Individuum
        // rangbasierte Selektion
        for (int i = 0; i < population.size(); i++) {

            Individuum currentIndividuum = population.get(i);
            currentIndividuum.setWahrscheinlichkeit(
                    berechneWahrscheinlichkeit(i + 1, population.size()));
            if (i > 0) {

                currentIndividuum
                        .setWahrsch_von(population.get(i - 1).getWahrsch_bis());
                currentIndividuum
                        .setWahrsch_bis(currentIndividuum.getWahrsch_von()
                                + currentIndividuum.getWahrscheinlichkeit());
            } else {
                // fuer das erste Element
                currentIndividuum.setWahrsch_von(0.0);
                currentIndividuum.setWahrsch_bis(
                        currentIndividuum.getWahrscheinlichkeit());
            }
        }
    }

    public static void bewertePopulation(List<Individuum> population) {
        for (Individuum individuum : population) {
            // Bei Nullstellen-Problem m�sste man jetzt zuerst eine L�sung
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

        // w�hle zuf�llig, welches Gen mutiert werden soll
        // wirft eine Zufallszahl zwischen 0 inkl. und genom.length exklusive
        int zufallGen = r.nextInt(genom.length);

        // Erzeuge eine Zufallszahl im definierten Wertebereich
        double zufallsZahl = (double) r.nextInt(
                Evolution.WERTEBEREICH_BIS - Evolution.WERTEBEREICH_VON + 1)
                + Evolution.WERTEBEREICH_VON;

        // mutiere das Gen
        genom[zufallGen] = genom[zufallGen] + zufallsZahl;
        // falls der Wertebereich �berschritten wird, sorge daf�r dass die
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
    public static Individuum rekombiniereElternIntermediaer(
            Individuum elternteil_A, Individuum elternteil_B) {

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
     * 
     * @param elternteil_A
     * @param elternteil_B
     * @return
     */
    public static Individuum rekombiniereElternArithmetisch(
            Individuum elternteil_A, Individuum elternteil_B) {

        double[] genom_A = elternteil_A.getGenom();
        double[] genom_B = elternteil_B.getGenom();

        double[] genom_Kind = new double[genom_A.length];

        // rekombiniere alle einzelnen Allele (zuf�lliger Wert zwischen den
        // beiden Werten)
        Random r = new Random();

        for (int i = 0; i < genom_A.length; i++) {

            double wertBis = Math.max(genom_A[i], genom_B[i]);
            double wertVon = Math.min(genom_A[i], genom_B[i]);
            int bound = 500;

            // next Double wirft Werte von 0 bis 1.0, wie strecken wir diesen
            // Wertebreich auf einen bestimmten Wertebereich?

            if (wertBis == wertVon) {
                genom_Kind[i] = wertVon;
            } else {
                double zufallszahl = r.nextDouble();
                genom_Kind[i] = zufallszahl * (wertBis - wertVon) + wertVon;
            }
        }

        return new Individuum(genom_Kind);
    }

    /**
     * Rouletteverfahren W�hlt zuf�llig ein Individuum aus einer Population
     * anhand des Rouletteverfahrens
     * 
     * @param population
     *            Die Population, aus denen ein Individuum gew�hlt werden soll
     * @return gew�hltes Individuum
     * @throws Exception
     */
    public static Individuum selektiereIndividuumRoulette(
            List<Individuum> population) throws Exception {

        // Erzeugen einer Zufallszahl zwischen 0 und 1.0
        // => muss so kompliziert gemacht werden, weil bei random die 1.0
        // exclusive ist
        Random r = new Random();
        double z = (double) r.nextInt(1000000001) / 1000000000;

        for (Individuum currentIndividuum : population) {
            // Roulette Auswahlverfahren
            if (currentIndividuum.getWahrsch_von() <= z
                    && z < currentIndividuum.getWahrsch_bis()) {
                currentIndividuum.setZufallszahlRoulette(z);
                return currentIndividuum;
                // currentIndividuum ist als Elternteil gew�hlt
            }
        }

        throw new Exception(
                "Es wurde kein Individuum mit der gew�rfelten Zufallszahl gefunden. Dies kann nicht sein.");
    }

    /**
     * Gibt rein zuf�llig ein Individuum aus der Population zur�ck.
     * 
     * @param population
     * @return
     */
    public static Individuum selektiereIndividuumZufaelligGleichverteilt(
            List<Individuum> population) {

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

    public static void printPopulationMitWahrscheinlichkeit(List<Individuum> population) {
        for (Individuum individuum : population) {
           printAusfuehrlichesIndividuum(individuum);
        }
    }
    /**
     * Gibt das Individuum auf der Konsole aus.
     * 
     * @param individuum
     */
    public static void printAusfuehrlichesIndividuum(Individuum individuum) {
        String format = "%.2f";
        String ausgabeString = "Genom: [";

        for (int i = 0; i < individuum.getGenom().length; i++) {
            ausgabeString += String.format(Locale.ENGLISH, format,
                    individuum.getGenom()[i]) + " ";
        }
        ausgabeString += "]";
        ausgabeString += " Fitness: " + String.format(Locale.ENGLISH, format,
                individuum.getFitness());
        ausgabeString += " - Wahrscheinlichkeit.: "
                + String.format(Locale.ENGLISH, format,
                        individuum.getWahrscheinlichkeit())
                + ", ";
        ausgabeString += "wird gewaehlt im Zufallszahlbereich: [" + String
                .format(Locale.ENGLISH, format, individuum.getWahrsch_von())
                + ", ";
        ausgabeString += String.format(Locale.ENGLISH, format,
                individuum.getWahrsch_bis()) + "]";

        System.out.println(ausgabeString);
    }
    
   
    
    public static void printIndividuum(Individuum individuum) {
        String format = "%.16f";
        String ausgabeString = "Genom: [";

        for (int i = 0; i < individuum.getGenom().length; i++) {
            ausgabeString += String.format(Locale.GERMAN, format,
                    individuum.getGenom()[i]) + " ";
        }
        ausgabeString += "]";
        ausgabeString += "\nFitness: " + String.format(Locale.GERMAN, format,
                individuum.getFitness());
        
        System.out.println(ausgabeString);
    }

    public static void printGenom(Individuum individuum) {
        String format = "%.2f";
        String ausgabeString = "Genom: [";

        for (int i = 0; i < individuum.getGenom().length; i++) {
            ausgabeString += String.format(Locale.ENGLISH, format,
                    individuum.getGenom()[i]) + " ";
        }
        ausgabeString += "]";

        System.out.println(ausgabeString);
    }

    private static void printGenomMitWahrscheinlichkeit(
            Individuum individuum) {
        String format = "%.2f";
        String ausgabeString = "Genom: [";

        for (int i = 0; i < individuum.getGenom().length; i++) {
            ausgabeString += String.format(Locale.ENGLISH, format,
                    individuum.getGenom()[i]) + " ";
        }
        ausgabeString += "] - Roulette-Zufallszahl: " + individuum.getZufallszahlRoulette();

        System.out.println(ausgabeString);
        
    }
    /***
     * Berechnung der Wahrsch. f�r rangbasierte Selektion f�r
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
     * 
     * @param text
     */
    public static void printCurrentStep(String text) {
        // System.out.println(text);
    }

}
