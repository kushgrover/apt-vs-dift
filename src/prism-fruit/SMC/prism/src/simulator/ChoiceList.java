//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

package simulator;

import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.doubles.DoubleListIterator;
import parser.State;
import parser.VarList;
import parser.ast.Command;
import parser.ast.Update;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLangException;

import java.util.ArrayList;
import java.util.List;

public class ChoiceList implements Choice
{
	protected String action;
	// List of multiple targets and associated info
	// Size of list is stored implicitly in target.length // TODO: no
	// TODO: convert to arrays?
	protected List<List<Update>> updates;
	protected DoubleList probability;
	protected List<Command> command;

	public ChoiceList(int n)
	{
		updates = new ArrayList<>(n);
		probability = new DoubleArrayList(n);
		command = new ArrayList<>(n);
	}

	// Set methods

	public void setAction(String action)
	{
		this.action = action;
	}

	public void setProbability(double probability)
	{
		setProbability(0, probability);
	}

	public void addProbability(double probability)
	{
		this.probability.add(probability);
	}

	public void setProbability(int i, double probability)
	{
		if (i < 0 || i >= size())
			return;
		this.probability.set(i, probability);
	}

	public void setCommand(Command command)
	{
		setCommand(0, command);
	}

	public void addCommand(Command command)
	{
		this.command.add(command);
	}

	public void setCommand(int i, Command command)
	{
		if (i < 0 || i >= size())
			return;
		this.command.set(i, command);
	}

	@Override
	public void scaleProbabilitiesBy(double d)
	{
		int i, n;
		n = size();
		for (i = 0; i < n; i++) {
			probability.set(i, probability.getDouble(i) * d);
		}
	}

	// Get methods

	@Override public int getModuleOrActionIndex()
	{
		return 0; // TODO
	}

	@Override public String getModuleOrAction()
	{
		return null; // TODO
	}

	public String getAction()
	{
		return action;
	}

	@Override public int size()
	{
		return probability.size();
	}

	@Override public String getUpdateString(int i, State currentState) throws PrismLangException
	{
		StringBuilder s = new StringBuilder("(");
		for (Update up : updates.get(i))
			s.append(up);
		s.append(")");
		return s.toString();
	}

	@Override public String getUpdateStringFull(int i)
	{
		return null;
	}

	public State computeTarget(State oldState) throws PrismLangException
	{
		return computeTarget(0, oldState);
	}

	public void computeTarget(State oldState, State newState) throws PrismLangException
	{
		computeTarget(0, oldState, newState);
	}

	@Override public State computeTarget(int i, State oldState) throws PrismLangException
	{
		if (i < 0 || i >= size())
			throw new PrismLangException("Choice does not have an element " + i);
		State newState = new State(oldState);
		for (Update up : updates.get(i))
			up.update(oldState, newState);
		return newState;
	}

	@Override public void computeTarget(int i, State oldState, State newState) throws PrismLangException
	{
		if (i < 0 || i >= size())
			throw new PrismLangException("Choice does not have an element " + i);
		for (Update up : updates.get(i))
			up.update(oldState, newState);
	}

	public double getProbability()
	{
		return getProbability(0);
	}

	@Override public double getProbability(int i)
	{
		if (i < 0 || i >= size())
			return -1;
			//throw new PrismLangException("Invalid grouped transition index " + i);
		return probability.getDouble(i);
	}

	@Override public double getProbabilitySum()
	{
		double sum = 0.0;
		DoubleListIterator iterator = probability.iterator();
		while (iterator.hasNext()) {
			sum += iterator.nextDouble();
		}
		return sum;
	}

	public Command getCommand()
	{
		return getCommand(0);
	}

	public Command getCommand(int i)
	{
		if (i < 0 || i >= size())
			return null;
			//throw new PrismLangException("Invalid grouped transition index " + i);
		return command.get(i);
	}

	@Override public int getIndexByProbabilitySum(double x)
	{
		int n = size();
		double d = 0.0;
		for (int i = 0; i < n; i++) {
			d += probability.getDouble(i);
			if (x < d) {
				return i - 1;
			}
		}
		return n;
	}

	@Override
	public void checkValid(ModelType modelType) throws PrismException
	{
		// TODO
	}

	@Override
	public void checkForErrors(State currentState, VarList varList) throws PrismException
	{
		// TODO
	}

	public String toString()
	{
		int i, n;
		boolean first = true;
		StringBuilder s = new StringBuilder();
		n = size();
		for (i = 0; i < n; i++) {
			if (first)
				first = false;
			else
				s.append(" + ");
			s.append(getProbability(i)).append(":").append(updates.get(i));
		}
		return s.toString();
	}
}
