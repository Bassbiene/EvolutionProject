package alg;

import java.util.Comparator;

/***
 * Ein Individuum mit kleinerer Fitness-Zahl kommt vor einem Individuum
 * mit großer Fitness-Zahl
 * @author Stefan
 *
 */
public class FitnessComparator implements Comparator<Individuum>{

	public int compare(Individuum i1, Individuum i2) {
				
		if (i1.getFitness() == null && i2.getFitness() == null) {
		      return 0;
		    }
		    if (i1.getFitness() == null) {
		      return 1;
		    }
		    if (i2.getFitness() == null) {
		      return -1;
		    }
		    return i1.getFitness().compareTo(i2.getFitness());
	}
}
