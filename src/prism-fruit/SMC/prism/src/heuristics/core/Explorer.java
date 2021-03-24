package heuristics.core;

import com.google.common.collect.Iterables;
import de.tum.in.naturals.set.NatBitSet;
import explicit.Distribution;
import it.unimi.dsi.fastutil.ints.IntIterable;
import prism.PrismException;

import java.util.List;

public interface Explorer
{
	void exploreState(int stateNumber) throws PrismException;

	boolean isStateExplored(int number);

	default int getInitialState()
	{
		return Iterables.getOnlyElement(getInitialStates());
	}

	IntIterable getInitialStates();

	int exploredStateCount();

	NatBitSet getExploredStates();

	List<Distribution> getChoices(int state);

	Distribution getChoice(int state, int action);
}
