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

package heuristics.nextstate;

import heuristics.nextstate.HeuristicNextStateFactory.NextState;
import heuristics.update.StateUpdate;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import parser.State;
import prism.ModelGenerator;
import prism.PrismException;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public abstract class HeuristicNextState
{
	static final Random random = new Random();

	protected ModelGenerator modelGenerator;
	protected StateUpdate stateUpdate;
	private State previousState;
	private int previousChoice;
	private int previousTransition;
	protected Map<State, Int2IntMap> stateToLastChosenTrans = new HashMap<>();

	public HeuristicNextState(ModelGenerator modelGenerator, StateUpdate stateUpdate)
	{
		this.modelGenerator = modelGenerator;
		this.stateUpdate = stateUpdate;
	}

	public abstract State sample(State state, int action, int depth) throws PrismException;

	public abstract State sample(State state, int action) throws PrismException;

	public abstract NextState getType();

	public State getPreviousState() throws PrismException
	{
		return previousState;
	}

	public int getPreviousChoice() throws PrismException
	{
		return previousChoice;
	}

	public int getPreviousTransitionOffset() throws PrismException
	{
		return previousTransition;
	}

	protected void updateChoice(State state, int action, int transitionOffset)
	{
		this.previousState = state;
		this.previousChoice = action;
		this.previousTransition = transitionOffset;
	}

	/**
	 * Sample according to the distribution
	 * @param state
	 * @param action
	 * @param dist
	 * @return
	 * @throws PrismException
	 */
	protected State sampleFromDist(State state, int action, double[] dist) throws PrismException
	{
		// Compute the cumulative sum (times 100)
		// eg. (0.2, 0.45, 0.35) -> (20, 65, 100)
		for (int i = 0; i < dist.length; i++) {
			if (i == 0) {
				dist[i] = Math.ceil(dist[i] * 100d);
			} else {
				dist[i] = Math.ceil(dist[i - 1] + dist[i] * 100d);
			}
		}

		int max = (int) Math.ceil(dist[dist.length - 1]);
		// max is 100 + (0..4) now

		if (max == 0) {
			final int offset = random.nextInt(dist.length);
			updateChoice(state, action, offset);
			return modelGenerator.computeTransitionTarget(action, offset);
		}
		int rand = random.nextInt(max) + 1;
		for (int i = 0; i < dist.length; i++) {
			if ((double) rand <= dist[i]) {
				updateChoice(state, action, i);
				return modelGenerator.computeTransitionTarget(action, i);
			}
		}
		throw new PrismException("Not sampling any state, this should never happen");
	}
}
