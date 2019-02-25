'use strict';

const {ipcRenderer, remote} = require('electron');
const storage = require('electron-json-storage');
const dialog = remote.require('electron').dialog;

// Request to go back to scenario list
const BacktoScenariosListButton =
    document.getElementById(('BacktoScenariosListButton'));
BacktoScenariosListButton.addEventListener(
    'click', function (event) {
        ipcRenderer.send("User-Requests-Index-View");
    }
);

// Get the database file
const dbFilePathSettingsButton =
    document.getElementById(('dbFilePathSettingsButton'));
dbFilePathSettingsButton.addEventListener(
    'click',
    function (event) {
    // Write where
    let dbFilePath = dialog.showOpenDialog({
        properties: ['openFile']
    });
    // set the database path
    storage.set(
        'dbFilePath',
        { 'dbFilePath': dbFilePath },
        function(error) {if (error) throw error;}
        );
    // re-load settings view
    ipcRenderer.send('User-Requests-Settings-View');
    }
);

storage.get(
    'dbFilePath',
    function(error, data) {
        if (error) throw error;
        console.log(data);
        document.getElementById('currentDBPath').innerHTML = data['dbFilePath'];
        }
    );


// Get the scenarios directory
const scenariosDirectorySettingsButton =
    document.getElementById(('scenariosDirectorySettingsButton'));
scenariosDirectorySettingsButton.addEventListener(
    'click',
    function (event) {
    // Write where
    let dbFilePath = dialog.showOpenDialog({
        properties: ['openDirectory']
    });
    // set the database path
    storage.set(
        'scenariosDirectory',
        { 'scenariosDirectory': dbFilePath },
        function(error) {if (error) throw error;}
        );
    // re-load settings view
    ipcRenderer.send('User-Requests-Settings-View');
    }
);

storage.get(
    'scenariosDirectory',
    function(error, data) {
        if (error) throw error;
        console.log(data);
        document.getElementById('currentScenariosDirectory').innerHTML = data['scenariosDirectory'];
        }
    );