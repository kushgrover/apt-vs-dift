//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Aistis Simaitis <aistis.aimaitis@cs.ox.ac.uk> (University of Oxford)
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

package strat;

import explicit.Model;
import prism.PrismLog;

import static com.google.common.base.Preconditions.checkArgument;

/**
 * Class to store a memoryless deterministic (MD) strategy, as a (Java) array of choice indices.
 */
public class MDStrategyArray extends MDStrategy
{
	// TODO Incorporate Strategy.Choice (isChoiceUndefined, isChoiceUnreachable, etc.)

	// Model associated with the strategy
	private explicit.NondetModel model;
	// Index of choice taken in each state (wrt model above)
	// Other possible values: -1 (unknown), -2 (arbitrary), -3 (unreachable)
	private int choices[];

	/**
	 * Creates an MDStrategyArray from an integer array of choices.
	 * The array may later be modified/delete - take a copy if you want to keep it.
	 */
	public MDStrategyArray(explicit.NondetModel model, int choices[])
	{
		checkArgument(choices.length == model.getNumStates());
		this.model = model;
		this.choices = choices;
	}

	// Methods for MDStrategy

	@Override
	public int getNumStates()
	{
		return model.getNumStates();
	}

	@Override
	public boolean isChoiceDefined(int s)
	{
		return choices[s] >= 0;
	}

	@Override
	public void setChoiceIndex(int state, int index) {
		choices[state] = index;
	}

	@Override
	public Strategy.Choice getChoice(int s)
	{
		switch (choices[s]) {
		case -1:
			return Choice.UNKNOWN;
		case -2:
			return Choice.ARBITRARY;
		case -3:
			return Choice.UNREACHABLE;
		default:
			return Choice.INDEX;
		}
	}

	@Override
	public int getChoiceIndex(int s)
	{
		return choices[s];
	}

	@Override
	public Object getChoiceAction(int s)
	{
		int c = choices[s];
		if (c >= 0) {
			return model.getAction(s, c);
		}
		if (c == -1) {
			return "?";
		}
		if (c == -2) {
			return "*";
		}
		return "-";
	}

	// Methods for Strategy

	@Override
	public void exportInducedModel(PrismLog out)
	{
		Model inducedModel = model.constructInducedModel(this);
		inducedModel.exportToPrismExplicitTra(out);
	}

	@Override
	public void exportDotFile(PrismLog out)
	{
		model.exportToDotFileWithStrat(out, null, this);
	}

	@Override
	public void clear()
	{
		choices = null;
	}
}
