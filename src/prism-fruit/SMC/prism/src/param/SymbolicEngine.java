//==============================================================================
//
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Hongyang Qu <hongyang.qu@cc.ox.ac.uk> (University of Oxford)
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

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.IntIterator;
import parser.State;
import parser.ast.Command;
import parser.ast.Expression;
import parser.ast.Module;
import parser.ast.ModulesFile;
import parser.ast.Update;
import parser.ast.Updates;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLangException;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class SymbolicEngine
{
	// Info on model being explored
	protected ModulesFile modulesFile;
	protected ModelType modelType;
	protected int numModules;
	// Synchronising action info
	protected List<String> synchs;
	protected int numSynchs;
	protected int synchModuleCounts[];

	// Temporary storage:

	// Element i,j of updateLists is a list of the updates from module i labelled with action j
	// (where j=0 denotes independent, otherwise 1-indexed action label)
	protected List<List<List<Updates>>> updateLists;
	// Bit j of enabledSynchs is set iff action j is currently enabled
	// (where j=0 denotes independent, otherwise 1-indexed action label)
	protected NatBitSet enabledSynchs;
	// Element j of enabledModules is a NatBitSet showing modules which enable action j
	// (where j=0 denotes independent, otherwise 1-indexed action label)
	protected NatBitSet[] enabledModules;
	protected ModelBuilder modelBuilder;
	protected FunctionFactory functionFactory;

	//  flag that suppresses warnings during calculateTransitions
	private boolean noWarnings = false;

	public SymbolicEngine(ModulesFile modulesFile, ModelBuilder modelBuilder, FunctionFactory functionFactory)
	{
		this.modelBuilder = modelBuilder;
		this.functionFactory = functionFactory;

		// Get info from model
		this.modulesFile = modulesFile;
		modelType = modulesFile.getModelType();
		numModules = modulesFile.getNumModules();
		synchs = modulesFile.getSynchs();
		numSynchs = synchs.size();

		// Compute count of number of modules using each synch action
		// First, compute and cache the synch actions for each of the modules
		List<HashSet<String>> synchsPerModule = new ArrayList<>(numModules);
		for (int i = 0; i < numModules; i++) {
			synchsPerModule.add(new HashSet<>(modulesFile.getModule(i).getAllSynchs()));
		}
		// Second, do the counting
		synchModuleCounts = new int[numSynchs];
		for (int j = 0; j < numSynchs; j++) {
			synchModuleCounts[j] = 0;
			String s = synchs.get(j);
			for (int i = 0; i < numModules; i++) {
				if (synchsPerModule.get(i).contains(s))
					synchModuleCounts[j]++;
			}
		}

		// Build lists/bitsets for later use
		updateLists = new ArrayList<>(numModules);
		for (int i = 0; i < numModules; i++) {
			updateLists.add(new ArrayList<>(numSynchs + 1));
			for (int j = 0; j < numSynchs + 1; j++) {
				updateLists.get(i).add(new ArrayList<>());
			}
		}
		enabledSynchs = NatBitSets.boundedSet(numSynchs + 1);
		enabledModules = new NatBitSet[numSynchs + 1];
		for (int j = 0; j < numSynchs + 1; j++) {
			enabledModules[j] = NatBitSets.boundedSet(numModules);
		}
	}

	public Expression getProbabilityInState(Updates ups, int i, State state) throws PrismLangException
	{
		Expression p = ups.getProbability(i);
		return (p == null) ? Expression.Double(1.0) : p;
	}

	static boolean hasMoreThanOneVariable(Expression exp)
	{
		int varNum = 0;
		try {
			varNum = exp.getAllVars().size();
			//System.out.println("varNum = " + varNum);
		} catch(PrismLangException e) {

		}

		if (varNum >1) {
			return true;
		} else {
			return false;
		}
	}

	public TransitionList calculateTransitions(State state, boolean noWarnings) throws PrismException
	{
		List<ChoiceListFlexi> chs;
		int i, j, k, l, n, count;
		TransitionList transitionList = new TransitionList();

		this.noWarnings = noWarnings;

		// Clear lists/bitsets
		transitionList.clear();
		for (i = 0; i < numModules; i++) {
			for (j = 0; j < numSynchs + 1; j++) {
				updateLists.get(i).get(j).clear();
			}
		}
		enabledSynchs.clear();
		for (i = 0; i < numSynchs + 1; i++) {
			enabledModules[i].clear();
		}

		// Calculate the available updates for each module/action
		// (update information in updateLists, enabledSynchs and enabledModules)
		for (i = 0; i < numModules; i++) {
			calculateUpdatesForModule(i, state);
		}
		//System.out.println("updateLists: " + updateLists);

		// Add independent transitions for each (enabled) module to list
		IntIterator iterator = enabledModules[0].iterator();
		while (iterator.hasNext()) {
			i = iterator.nextInt();
			for (Updates ups : updateLists.get(i).get(0)) {
				transitionList.add(processUpdatesAndCreateNewChoice(-(i + 1), ups, state));
			}
		}
		// Add synchronous transitions to list
		chs = new ArrayList<>();
		iterator = enabledSynchs.iterator();
		while (iterator.hasNext()) {
			i = iterator.nextInt();
			if (i == 0) {
				continue;
			}
			chs.clear();
			// Check counts to see if this action is blocked by some module
			if (enabledModules[i].size() < synchModuleCounts[i - 1])
				continue;
			// If not, proceed...
			IntIterator moduleIterator = enabledModules[i].iterator();
			while (moduleIterator.hasNext()) {
				j = moduleIterator.nextInt();
				count = updateLists.get(j).get(i).size();
				// Case where there is only 1 Updates for this module
				if (count == 1) {
					Updates ups = updateLists.get(j).get(i).get(0);
					// Case where this is the first Choice created
					if (chs.isEmpty()) {
						chs.add(processUpdatesAndCreateNewChoice(i, ups, state));
					}
					// Case where there are existing Choices
					else {
						// Product with all existing choices
						for (ChoiceListFlexi ch : chs) {
							processUpdatesAndAddToProduct(ups, state, ch);
						}
					}
				}
				// Case where there are multiple Updates (i.e. local nondeterminism)
				else {
					// Case where there are no existing choices
					if (chs.isEmpty()) {
						for (Updates ups : updateLists.get(j).get(i)) {
							chs.add(processUpdatesAndCreateNewChoice(i, ups, state));
						}
					}
					// Case where there are existing Choices
					else {
						// Duplicate (count-1 copies of) current Choice list
						n = chs.size();
						for (k = 0; k < count - 1; k++)
							for (l = 0; l < n; l++)
								chs.add(new ChoiceListFlexi(chs.get(l)));
						// Products with existing choices
						for (k = 0; k < count; k++) {
							Updates ups = updateLists.get(j).get(i).get(k);
							for (l = 0; l < n; l++) {
								processUpdatesAndAddToProduct(ups, state, chs.get(k * n + l));
							}
						}
					}
				}
			}
			// Add all new choices to transition list
			for (ChoiceListFlexi ch : chs) {
				transitionList.add(ch);
			}
		}

		// Check validity of the computed transitions
		// (not needed currently)
		//transitionList.checkValid(modelType);

		// Check for errors (e.g. overflows) in the computed transitions
		//transitionList.checkForErrors(state, varList);

		//System.out.println(transitionList);
		return transitionList;
	}

	// Private helpers

	/**
	 * Determine the enabled updates for the 'm'th module from (global) state 'state'.
	 * Update information in updateLists, enabledSynchs and enabledModules.
	 * @param m The module index
	 * @param state State from which to explore
	 */
	protected void calculateUpdatesForModule(int m, State state) throws PrismLangException
	{
		Module module;
		Command command;
		int i, j, n;

		module = modulesFile.getModule(m);
		n = module.getNumCommands();
		for (i = 0; i < n; i++) {
			command = module.getCommand(i);
			if (command.getGuard().evaluateBoolean(state)) {
				j = command.getSynchIndex();
				updateLists.get(m).get(j).add(command.getUpdates());
				enabledSynchs.set(j);
				enabledModules[j].set(m);
			}
		}
	}

	/**
	 * Create a new Choice object (currently ChoiceListFlexi) based on an Updates object
	 * and a (global) state. Check for negative probabilities/rates.
	 * @param moduleOrActionIndex Module/action for the choice, encoded as an integer (see Choice)
	 * @param ups The Updates object
	 * @param state Global state
	 */
	private ChoiceListFlexi processUpdatesAndCreateNewChoice(int moduleOrActionIndex, Updates ups, State state) throws PrismLangException
	{
		ChoiceListFlexi ch;
		List<Update> list;
		int i, n;
		Expression p;

		// Create choice and add all info
		ch = new ChoiceListFlexi();
		ch.setModuleOrActionIndex(moduleOrActionIndex);
		n = ups.getNumUpdates();
		for (i = 0; i < n; i++) {
			// Compute probability/rate
			p = getProbabilityInState(ups, i, state);
			int[] varMap = new int[state.varValues.length];
			for (int var = 0; var < varMap.length; var++) {
				varMap[var] = var;
			}
			p = (Expression) p.deepCopy().evaluatePartially(state, varMap);
			list = new ArrayList<>();
			list.add(ups.getUpdate(i));

			try {
				Function pFn = modelBuilder.expr2function(functionFactory, p);
				if (pFn.isZero()) {
					// function for probability / rate is zero, don't add the corresponding transition
					if (!noWarnings)
						modelBuilder.getLog().printWarning("Update has zero " + (modelType.continuousTime() ? "rate" : "probability") + " (" + p + (p.hasPosition() ? ", " + p.getBeginString() : "") +")");
					continue;
				}
				ch.add(pFn, list);
			} catch (PrismException e) {
				throw new PrismLangException(e.getMessage());
			}
		}

		return ch;
	}

	/**
	 * Create a new Choice object (currently ChoiceListFlexi) based on the product
	 * of an existing ChoiceListFlexi and an Updates object, for some (global) state.
	 * If appropriate, check probabilities sum to 1 too.
	 * @param ups The Updates object
	 * @param state Global state
	 * @param ch The existing Choices object
	 */
	private void processUpdatesAndAddToProduct(Updates ups, State state, ChoiceListFlexi ch) throws PrismLangException
	{
		// Create new choice (action index is 0 - not needed)
		ChoiceListFlexi chNew = processUpdatesAndCreateNewChoice(0, ups, state);
		// Build product with existing
		ch.productWith(chNew);
	}
}
