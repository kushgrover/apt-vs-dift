//==============================================================================
//
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	* Mateusz Ujma <mateusz.ujma@cs.ox.ac.uk> (University of Oxford)
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

package heuristics.update;

import heuristics.CachedModelGenerator;
import heuristics.init.ValueInit;
import heuristics.search.StateValue;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.objects.Object2BooleanArrayMap;
import it.unimi.dsi.fastutil.objects.Object2BooleanMap;
import parser.State;
import parser.Values;
import parser.ast.Expression;
import prism.ModelGenerator;
import prism.PrismException;
import prism.PrismUtils;

public abstract class StateUpdate
{
	private final ValueInit valueInit = new ValueInit(this);
	private final Object2BooleanMap<State> cacheStateToTarget = new Object2BooleanArrayMap<>();
	private final Object2BooleanMap<State> cacheStateToSelfLoop = new Object2BooleanArrayMap<>();
	private final Object2BooleanMap<State> cacheStateToZero = new Object2BooleanArrayMap<>();
	protected final CachedModelGenerator modelGenerator;
	protected final int bound;
	protected final boolean min;
	protected final double epsilon;
	private Values constantValues;
	private Expression target;
	private State initialState;
	private final StateValueContainer stateValueContainer;
	private IntSet coalition;

	public StateUpdate(CachedModelGenerator modelGenerator, StateValueContainer container, int bound, boolean min, double epsilon) throws PrismException
	{
		this.modelGenerator = modelGenerator;
		this.bound = bound;
		this.stateValueContainer = container;
		if (modelGenerator != null) {
			this.initialState = this.modelGenerator.getInitialState();
		}
		this.min = min;
		this.epsilon = epsilon;
	}

	public abstract int update(State state, int depth) throws PrismException;

	public abstract int update(State state) throws PrismException;

	public boolean isTarget(State state) throws PrismException
	{
		// In case target doesn't exist
		if (target == null) {
			return false;
		}
		if (cacheStateToTarget.containsKey(state)) {
			return cacheStateToTarget.getBoolean(state);
		}
		modelGenerator.exploreState(state);
		Expression evTarget = target.deepCopy();
		boolean cachedValue = evTarget.evaluateBoolean(constantValues, state);
		cacheStateToTarget.put(state, cachedValue);
		return cachedValue;
	}

	public boolean isSelfLoop(State state) throws PrismException
	{
		if (cacheStateToSelfLoop.containsKey(state)) {
			return cacheStateToSelfLoop.getBoolean(state);
		}
		modelGenerator.exploreState(state);
		int numChoices = modelGenerator.getNumChoices();
		for (int i = 0; i < numChoices; i++) {
			int trans = modelGenerator.getNumTransitions(i);
			for (int j = 0; j < trans; j++) {
				State successor = modelGenerator.computeTransitionTarget(i, j);
				if (!successor.equals(state)) {
					cacheStateToSelfLoop.put(state, false);
					return false;
				}
			}
		}
		cacheStateToSelfLoop.put(state, true);
		return true;
	}

	public boolean isZero(State state) {
		return cacheStateToZero.getOrDefault(state, false);
	}

	public boolean isZeroState(State state) throws PrismException
	{
		return cacheStateToZero.containsKey(state);
	}

	/**
	 * Returns true if setTarget(state) has been executed previously
	 */
	public boolean isTargetState(State state) throws PrismException
	{
		return cacheStateToTarget.getOrDefault(state, false);
	}

	public boolean visited(State state)
	{
		return getQValue(state) == null;
	}

	public void setTarget(State state, boolean t) {
		setTarget(state, t, -1);
	}

	public void setTarget(State state, boolean t, int depth)
	{
		cacheStateToTarget.put(state, t);
		if (t) {
			setQValue(state, new StateValue(1d, 1d), depth);
		}
	}

	public void setZero(State state, boolean t)
	{
		setZero(state, t, -1);
	}

	public void setZero(State state, boolean t, int depth)
	{
		cacheStateToZero.put(state, t);
		if (t) {
			setQValue(state, new StateValue(0d, 0d), depth);
		}
	}

	public int getBound()
	{
		return bound;
	}

	public void setQValue(State state, StateValue val, int depth)
	{
		stateValueContainer.setQValue(state, val, depth);
	}

	public StateValue getQValue(State state, int depth)
	{
		return stateValueContainer.getQValue(state, depth);
	}

	public void setQValue(State state, StateValue val)
	{
		stateValueContainer.setQValue(state, val);
	}

	public StateValue getQValue(State state)
	{
		return stateValueContainer.getQValue(state);
	}

	public int getCurrentAction(State state)
	{
		Integer cA = stateValueContainer.getCurrentAction(state);
		if (cA == null) {
			return 0;
		}
		return cA;
	}

	public void setCurrentAction(State state, int a)
	{
		stateValueContainer.setCurrentAction(state, a);
	}

	public void setTargetOrBoundValue(State state, int depth) throws PrismException
	{
		if (isTarget(state)) {
			setQValue(state, new StateValue(1.0, 1.0), depth);
			return;
		}
		if (depth >= bound) {
			setQValue(state, new StateValue(0.0, 0.0), depth);
			return;
		}
	}

	public void setTargetOrZeroValue(State state) throws PrismException
	{
		if (isTarget(state)) {
			setQValue(state, new StateValue(1.0, 1.0));
		} else if (isZeroState(state)) {
			setQValue(state, new StateValue(0.0, 0.0));
		}
	}

	public int getStateValueNumber()
	{
		int stateValuesNumber = 0;
		/*Iterator<State> it = qValues.keySet().iterator();
		while(it.hasNext()) {
			State s = it.next();
			stateValuesNumber += qValues.get(s).size();
			if(qValues.get(s).size() == 0) {
				System.out.println("s " + s);
				System.exit(1);
			}
		}*/
		return stateValuesNumber;
	}

	public int getNumExplored()
	{
		return stateValueContainer.getSize();
	}

	public int getNumAllValues()
	{
		return stateValueContainer.getNumAllValues();
	}

	public boolean same(double d1, double d2)
	{
		return PrismUtils.doublesAreClose(d1, d2, epsilon, true);
	}

	public void setConstantValues(Values cV)
	{
		constantValues = cV;
	}

	public void setTarget(Expression e)
	{
		this.target = e;
	}

	public Expression getTarget()
	{
		return this.target;
	}

	public void setCoalition(IntSet coalition)
	{
		this.coalition = coalition;
	}

	public IntSet getCoalition()
	{
		return this.coalition;
	}

	/* TODO: PORT
	public int getPlayer(State s) throws PrismException
	{
		modelGenerator.exploreState(s);
		return coalition.contains(modelGenerator.getPlayerForState()) ? STPGExplicit.PLAYER_1 : STPGExplicit.PLAYER_2;
	}
	*/

	public ModelGenerator getModelExplorer()
	{
		return this.modelGenerator;
	}

	public StateValueContainer getStateValueContainer()
	{
		return this.stateValueContainer;
	}

	/**
	 * Returns the weighted average of the lower and upper bounds of the successors
	 */
	protected StateValue getActionValueBounds(State state, int a, int depth) throws PrismException
	{
		modelGenerator.exploreState(state);
		int numTransitions = modelGenerator.getNumTransitions(a);
		double lowerBound = 0d;
		double upperBound = 0d;
		for (int i = 0; i < numTransitions; i++) {
			modelGenerator.exploreState(state);
			final double prob = modelGenerator.getTransitionProbability(a, i);
			final State successor = modelGenerator.computeTransitionTarget(a, i);
			final StateValue successorStateValue;
			if (depth == -1) {
				successorStateValue = getQValue(successor);
			} else {
				successorStateValue = getQValue(successor, depth);
			}

			if (successorStateValue == null) {
				if (depth == -1) {
					lowerBound += prob * valueInit.getLowerBoundValue(modelGenerator, successor);
					upperBound += prob * valueInit.getUpperBoundValue(modelGenerator, successor);
				} else {
					lowerBound += prob * valueInit.getLowerBoundValue(modelGenerator, successor, depth);
					upperBound += prob * valueInit.getUpperBoundValue(modelGenerator, successor, depth);
				}
			} else {
				upperBound += prob * successorStateValue.upperBound;
				lowerBound += prob * successorStateValue.lowerBound;
			}
		}

		return new StateValue(lowerBound, upperBound);
	}

	protected StateValue getActionValueBounds(State state, int a) throws PrismException
	{
		return getActionValueBounds(state, a, -1);
	}
}
