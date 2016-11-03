package alg;

public class Individuum {

	public Individuum(double[] allele){
		this.genom = allele;
	}
	
	private Double fitness;	
	private double[] genom;
	private Double wahrscheinlichkeit;
	public Double getWahrsch_von() {
		return wahrsch_von;
	}

	public void setWahrsch_von(Double wahrsch_von) {
		this.wahrsch_von = wahrsch_von;
	}

	public Double getWahrsch_bis() {
		return wahrsch_bis;
	}

	public void setWahrsch_bis(Double wahrsch_bis) {
		this.wahrsch_bis = wahrsch_bis;
	}

	private Double wahrsch_von;
	private Double wahrsch_bis;
	
	public Double getWahrscheinlichkeit() {
		return wahrscheinlichkeit;
	}

	public void setWahrscheinlichkeit(Double wahrscheinlichkeit) {
		this.wahrscheinlichkeit = wahrscheinlichkeit;
	}

	public Double getFitness() {
		return fitness;
	}

	public void setFitness(Double fitness) {
		this.fitness = fitness;
	}

	public double[] getGenom() {
		return genom;
	}

	public void setGenom(double[] allele) {
		this.genom = allele;
	}


	
}
