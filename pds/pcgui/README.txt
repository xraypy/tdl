=== Notes on wxGUI use. ===

 * The wxGUI modules are run using wxPython and Python card.
   Both these packages need to be installed in the usual
   python site-packages directory for this to work.

 * Note these modules fail with wx2.8 (work fine with wx2.6).
   In the later versions hitting enter in a field acts like a tab
   in the old version - so moved focus from one widget to another.
   This is due to a change in the default behaviour of wx in the
   latter versions.  The fix is to edit the following in the PythonCard
   module to add wx.WANTS_CHAR and remove wx.TAB_TRAVERSAL from the style
   argument to wx.Panel.__init__() and wx.Frame.__init__():

   * model.py, line 671:
      'style | wx.NO_FULL_REPAINT_ON_RESIZE | wx.WANTS_CHARS,'
   * model.py, line 1345:
      'style=wx.WANTS_CHARS,'
   * widget.py, line 427:
      'style=wx.WANTS_CHARS | wx.NO_FULL_REPAINT_ON_RESIZE)'

 * The data used to render a window/set its properties are in
   the .rscr.py files. These define the various "components" that
   can then be handled in the wrapping class

 * Starting up the python card rescourceEditor and opening
   these files lets you sketch a window and add components,
   tags, values etc. using a wysiswyg (pretty nice)

 * Once launched the corresponding .py file handles events
   through a simple "on_component_event" syntax

=== Child windows ===

 * A child window can get data from the parent with a
   call to self.getParent().XXX call.  So the child is responsible
   for getting all the data.

 * Not sure about race conditions or blocking etc but so far this seems
  to work fine (eg the child can hold a pointer to the shell).   
 
=== Handling child windows: ===

 * Option 1:

    * currently the application windows are hard coded into the app_menu.menu_Apps class
    * this limts flexibilty / control at startup over what gets loaded

 * Option 2 (better?): 
    * Drop all the python card wx files into the wxGUI directory
      (as is done now).  

    * Create a startup file for wxGUI has name matches between
      tdl modules and the python card wxXXX.rscr.py and wxXXX.py

    * So at startup read the file, chomp through the list and
      dynamically create everything. e.g. this could be done in an init method in app_menu.menu_Apps

 {{{
   if tdl.module.isLoaded(xrf):
       execute code to add the xrf item to the app_menu.menu_Apps class
   end
 }}} 
