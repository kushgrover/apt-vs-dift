package simulator.sampler;

import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import java.util.HashSet;
import java.util.Set;

/**
 * Represents a candidate.
 *
 * @author Przemysław Daca
 */
public class Candidate<T>
{
	// path index where the candidate was born
	private final int birth;
	public int firstentry;
	protected Set<T> states;
	// how many times a state occurs after the birth, null means that all entries are zero
	private Object2IntMap<T> counter;

	// The candidate is known to not to be bottom
	private boolean incomplete = false;
	// The candidate is not a SCCC
	private boolean trivial = true;

	public Candidate(int idx, T st)
	{
		birth = idx;
		states = new HashSet<>();
		states.add(st);
		firstentry = idx;
	}

	/**
	 * Count occurrence of candidate's state.
	 *
	 * @param st
	 */
	public void countState(T st)
	{
		if (counter == null) {
			counter = new Object2IntOpenHashMap<>();
			counter.defaultReturnValue(0);
		}

		counter.put(st, counter.getInt(st) + 1);
	}

	public boolean isTrivial()
	{
		return trivial;
	}

	public void setNonTrivial()
	{
		trivial = false;
	}

	public Integer getCount(T o)
	{
		int indx = counter.getOrDefault(o, -1);
		if (indx == -1) {
			return null;
		}
		return indx;
	}

	public boolean contains(T st)
	{
		return states.contains(st);
	}

	public Set<T> getStates()
	{
		return states;
	}

	public int getBirthIdx()
	{
		return birth;
	}

	/**
	 * Merge with another candidate.
	 *
	 * @param old
	 */
	public void mergeWith(Candidate<T> old)
	{
		states.addAll(old.states);
		if (old.firstentry < firstentry) {
			firstentry = old.firstentry;
		}
	}

	public String toString()
	{
		return "Candidate birth=" + birth + ", states=" + states;
	}

	/**
	 * Returns true iff the candidate is k-strong.
	 *
	 * @param k
	 * @return
	 */
	public boolean isStrong(int k)
	{
		if (counter == null) {
			return false;
		}

		for (T st : states) {
			if (counter.getInt(st) < k) {
				return false;
			}
		}

		return true;
	}

	public void setIncomplete()
	{
		incomplete = true;
	}

	public boolean isIncomplete()
	{
		return incomplete;
	}

}