//==============================================================================
//
//	Copyright (c) 2002-
//	Authors:
//	* Andrew Hinton <ug60axh@cs.bham.ac.uk> (University of Birmingham)
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
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

package userinterface;

import userinterface.util.GUIEvent;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JMenu;
import javax.swing.KeyStroke;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

public class GUIFileMenu extends GUIPlugin
{
    private JMenu menu;
    private Action exitAction;

    public GUIFileMenu(GUIPrism pr)
    {
        super(pr, false);
        initComponents();
    }

	@Override public void takeCLArgs(String args[])
	{
	}

    @Override public boolean displaysTab()
    {
        return false;
    }

    @Override public javax.swing.JMenu getMenu()
    {
        return menu;
    }

    @Override public OptionsPanel getOptions()
    {
        return null;
    }

    @Override public String getTabText()
    {
        return "";
    }

    @Override public javax.swing.JToolBar getToolBar()
    {
        return null;
    }

    @Override public boolean processGUIEvent(GUIEvent e)
    {
		return false;
    }

    private void initComponents()
    {
        menu = new JMenu("File");
        menu.setMnemonic('F');
        exitAction = new AbstractAction()
        {
            @Override public void actionPerformed(ActionEvent e)
            {
                getGUI().exit();
            }
        };
        exitAction.putValue(Action.LONG_DESCRIPTION, "Exits the application");
        //exitAction.putValue(Action.SHORT_DESCRIPTION, "Exit");
        exitAction.putValue(Action.MNEMONIC_KEY, new Integer(KeyEvent.VK_X));
        exitAction.putValue(Action.NAME, "Exit");
        exitAction.putValue(Action.SMALL_ICON, GUIPrism.getIconFromImage("smallExit.png"));
        exitAction.putValue(Action.ACCELERATOR_KEY, KeyStroke.getKeyStroke(KeyEvent.VK_Q, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
        menu.add(exitAction);
    }

    @Override public void notifySettings(prism.PrismSettings settings)
    {}

}
