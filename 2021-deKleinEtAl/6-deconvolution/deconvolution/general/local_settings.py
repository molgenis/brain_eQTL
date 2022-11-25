"""
File:         local_settings.py
Created:      2020/03/19
Last Changed: 2020/03/20
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import json
import os

# Third party imports.

# Local application imports.


class LocalSettings:
    """
    LocalSettings: class that utilizes the default_settings.json.
    """
    class __LocalSettings:
        """
        Singleton class that enables the extraction of data from the
        default_settings.json file.
        """
        def __init__(self, project_root, settings_file):
            """
            Initializer method for the LocalSettings class.
            """
            self.settings = None
            self.settings_path = os.path.join(project_root,
                                              'settings',
                                              settings_file + ".json")
            self.load_settings()

        def load_settings(self):
            """
            Method to load the default_settings.json and
            saves it to self.settings.
            """
            settings_content = open(self.settings_path)
            settings = json.loads(settings_content.read())
            settings_content.close()
            self.settings = settings

        def is_settings_loaded(self):
            """
            Method to check if settings are loaded.

            :return Boolean: true of loaded and false when not.
            """
            return self.settings is not None

        def get_setting(self, setting_key):
            """
            Method to return the value of a local_setting key.

            :param setting_key: string, key in settings.
            :return: str/dict, value of setting_key. If the key does not
                     exist; return None.
            """
            if self.settings is not None and setting_key in self.settings:
                value = self.settings[setting_key]
            else:
                value = None

            return value

        def get_all_settings(self):
            """
            Method to return all the local settings file.

            :return self.settings: dict, dictionary containing all local
                    settings.
            """
            return self.settings

        def __str__(self):
            """
            __str__ method to return a repr of self.

            :return: repr(self)
            """
            return repr(self)

    instance = None

    def __new__(cls, project_root, settings_file):
        """
        Class method that handels object creation. It is responsible for
        returning the LocalSettings class instance. If the instance is None
        (i.e. not there), create a new one and safe / return that.

        :return: instance
        """
        if not LocalSettings.instance:
            LocalSettings.instance = LocalSettings.__LocalSettings(project_root,
                                                                   settings_file)
        return LocalSettings.instance

    def __init__(self, project_root, settings_file):
        """
        Initializer method setting the LocalSettings class instance.
        """
        if not LocalSettings.instance:
            LocalSettings.instance = LocalSettings.__LocalSettings(project_root,
                                                                   settings_file)

    def __getattr__(self, name):
        """
        Method that is called when an attribute lookup has not found the
        attribute in the usual places. This returns the class instance from
        the instance variable.

        :param name: str
        :return: str
        """
        return getattr(self.instance, name)
