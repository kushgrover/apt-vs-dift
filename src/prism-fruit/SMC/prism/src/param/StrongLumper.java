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

package param;

import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntListIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.function.IntConsumer;

/**
 * Strong bisimulation lumper.
 *
 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
 */
final class StrongLumper extends Lumper {

	/**
	 * Construct a new strong bisimulation lumper.
	 *
	 * @param origPmc Markov chain to construct lumper for
	 */
	StrongLumper(MutablePMC origPmc) {
		super(origPmc);
	}

	/**
	 * Construct the strong bisimulation signature of given block.
	 * The signature is a mapping of blocks to the probability to move
	 * from the given state to any state of the block.
	 *
	 * @param state state to compute signature of
	 * @return signature of this state
	 */
	private Map<IntSet, Function> stateSignature(int state)
	{
		Map<IntSet, Function> signature = new HashMap<>();
		IntListIterator toStateIter = origPmc.transitionTargets.get(state).listIterator();
		ListIterator<Function> toProbIter = origPmc.transitionProbs.get(state).listIterator();
		while (toStateIter.hasNext()) {
			int toState = toStateIter.nextInt();
			Function toStateProb = toProbIter.next();
			IntSet toBlock = partition.getStateBlock(toState);
			Function toBlockProb = signature.get(toBlock);
			if (toBlockProb == null) {
				toBlockProb = origPmc.getFunctionFactory().getZero();
			}
			toBlockProb = toBlockProb.add(toStateProb);
			signature.put(toBlock, toBlockProb);
		}
		return signature;
	}

	/**
	 * Refines a given block to a list of new blocks for strong bisimulation.
	 * Each block will consist of the states with the same signature.
	 *
	 * @param oldBlock block to refine
	 * @param newBlocks list of new blocks generated
	 */
	@Override
	protected void refineBlock(IntSet oldBlock, List<IntSet> newBlocks) {
		Map<Map<IntSet, Function>, IntSet> signatures = new HashMap<>();
		oldBlock.forEach((IntConsumer) state -> {
			Map<IntSet, Function> signature = stateSignature(state);
			IntSet newBlock = signatures.computeIfAbsent(signature, k -> new IntOpenHashSet());
			newBlock.add(state);
		});
		newBlocks.addAll(signatures.values());
	}

	/**
	 * Build the strong bisimulation quotient from the blocks computed.
	 * Transition probabilities are simply derived from the signature
	 * of an arbitrary state for each block.
	 */
	@Override
	protected void buildQuotient() {
		optPmc = new MutablePMC(origPmc.getFunctionFactory(), blocks.size(), origPmc.isUseRewards(), origPmc.isUseTime());
		for (int newState = 0; newState < blocks.size(); newState++) {
			IntSet fromBlock = blocks.get(newState);
			int someOriginalState = fromBlock.iterator().nextInt();
			Map<IntSet, Function> signature = stateSignature(someOriginalState);
			for (Entry<IntSet, Function> entry : signature.entrySet()) {
				optPmc.addTransition(newState, blockToNumber.getInt(entry.getKey()), entry.getValue());
			}
			if (origPmc.isUseRewards()) {
				optPmc.setReward(newState, origPmc.getReward(someOriginalState));
				if (optPmc.getReward(newState).equals(origPmc.getFunctionFactory().getZero())) {
					optPmc.setTargetState(newState, true);
				}
				if (origPmc.isUseTime()) {
					optPmc.setTime(newState, origPmc.getTime(someOriginalState));
				}
			} else {
				optPmc.setTargetState(newState, origPmc.isTargetState(someOriginalState));
			}
			IntIterator iterator = fromBlock.iterator();
			while (iterator.hasNext()) {
				int oldState = iterator.nextInt();
				if (origPmc.isInitState(oldState)) {
					optPmc.setInitState(newState, true);
				}
			}
		}
	}

}
