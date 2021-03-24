package heuristics.update;

import heuristics.search.StateValue;
import parser.State;

/**
 * Created by tlm on 20/01/17.
 */
public interface StateValueContainer
{
	default boolean setQValue(State state, StateValue value)
	{
		return setQValue(state, value, -1);
	}

	void setCurrentAction(State state, int action);

	int getCurrentAction(State state);

	boolean setQValue(State state, StateValue value, int depth);

	default StateValue getQValue(State state)
	{
		return getQValue(state, -1);
	}

	StateValue getQValue(State state, int depth);

	int getSize();

	int getNumAllValues();
}
