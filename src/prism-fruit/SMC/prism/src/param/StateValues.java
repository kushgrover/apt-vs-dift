//==============================================================================
//
//	Copyright (c) 2013-
//	Authors:
//	* Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
//
//------------------------------------------------------------------------------
//
//	This file is part of PRISM.
//
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//==============================================================================

// TODO allow null values in entries

package param;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;

import java.util.ArrayList;

/**
 * Class to assign a value to each state of a model.
 * Used by {@code RegionValues}.
 *
 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
 * @see RegionValues
 */
public final class StateValues
{
	/** assigns values to each state of the model */
	private final ArrayList<StateValue> values;
	/** initial state of the model */
	private final int initState;

	/**
	 * Constructs new set of state values.
	 * Each state is mapped to {@code null}.
	 *
	 * @param numStates number of states of model
	 * @param initState initial state of the model
	 */
	public StateValues(int numStates, int initState)
	{
		values = new ArrayList<>(numStates);
		for (int state = 0; state < numStates; state++) {
			values.add(state, null);
		}
		this.initState = initState;
	}

	/**
	 * Constructs new set of state values.
	 * Each state is mapped to the given value.
	 *
	 * @param numStates number of states of the model
	 * @param initState initial state of the model
	 * @param value value to map all states to
	 */
	public StateValues(int numStates, int initState, StateValue value)
	{
		this(numStates, initState);
		for (int state = 0; state < numStates; state++) {
			values.set(state, value);
		}
	}

	/**
	 * Constructs new set of state values.
	 * Each state is mapped to the given value.
	 *
	 * @param numStates number of states of the model
	 * @param initState initial state of the model
	 * @param value value to map all states to
	 */
	public StateValues(int numStates, int initState, boolean value)
	{
		this(numStates, initState, new StateBoolean(value));
	}

	@Override
	public String toString()
	{
		return values.get(initState).toString();
	}

	@Override
	public boolean equals(Object obj)
	{
		if (!(obj instanceof StateValues)) {
			return false;
		}

		StateValues result = (StateValues) obj;

		for (int i = 0; i < values.size(); i++) {
			if (!values.get(i).equals(result.values.get(i))) {
				return false;
			}
		}

		return true;
	}

	@Override
	public int hashCode()
	{
		int hash = 0;

		for (int i = 0; i < values.size(); i++) {
			hash = values.get(i).hashCode() + (hash << 6) + (hash << 16) - hash;
		}

		return hash;
	}

	/**
	 * Get value of given state.
	 *
	 * @param state state to get value of
	 * @return value of given state
	 */
	public StateValue getStateValue(int state)
	{
		return values.get(state);
	}

	/**
	 * Set value of given state.
	 *
	 * @param state state to set value of
	 * @param value value to set for state
	 */
	public void setStateValue(int state, StateValue value)
	{
		values.set(state, value);
	}

	/**
	 * Set value of given state.
	 *
	 * @param state state to set value of
	 * @param value value to set for state
	 */
	public void setStateValue(int state, boolean value)
	{
		values.set(state, new StateBoolean(value));
	}

	/**
	 * Get value of given state as rational function.
	 * If the value of the state is not a function, this will lead to an error.
	 *
	 * @param state state to get value of
	 * @return value of the state as a function
	 */
	public Function getStateValueAsFunction(int state)
	{
		return (Function) values.get(state);
	}

	/**
	 * Get value of given state as boolean.
	 * If the value of the state is not a boolean, this will lead to an error.
	 *
	 * @param state state to get value of
	 * @return value of the state as a boolean
	 */
	public boolean getStateValueAsBoolean(int state)
	{
		return ((StateBoolean) values.get(state)).getValue();
	}

	/**
	 * Get value of initial state as rational function.
	 * If the value of the initial state is not a function, this will lead to an error.
	 *
	 * @return value of the initial state as a function
	 */
	public Function getInitStateValueAsFunction()
	{
		return (Function) values.get(initState);
	}

	/**
	 * Get value of initial state as boolean.
	 * If the value of the initial state is not a boolean, this will lead to an error.
	 *
	 * @return value of the initial state as a boolean
	 */
	public boolean getInitStateValueAsBoolean()
	{
		return ((StateBoolean) values.get(initState)).getValue();
	}

	/**
	 * Returns number of states of the model.
	 *
	 * @return number of states of the model
	 */
	public int getNumStates()
	{
		return values.size();
	}

	/**
	 * Converts this state value assignment to a bitset.
	 * For this to work, all states must be mapped to booleans.
	 *
	 * @return bitset representing this state value assignment
	 */
	public NatBitSet toBitSet()
	{
		NatBitSet result = NatBitSets.boundedSet(values.size());
		for (int state = 0; state < values.size(); state++) {
			result.set(state, getStateValueAsBoolean(state));
		}
		return result;
	}

	/**
	 * Instantiates the value to which each state is mapped at the given point.
	 * For this to work, all states must be mapped to rational functions.
	 *
	 * @param point point to instantiate state values
	 * @return array of {@code BigRational}s mapping each state to evaluated value
	 */
	public BigRational[] instantiate(Point point)
	{
		BigRational[] result = new BigRational[values.size()];
		for (int state = 0; state < values.size(); state++) {
			result[state] = this.getStateValueAsFunction(state).evaluate(point);
		}
		return result;
	}
}
