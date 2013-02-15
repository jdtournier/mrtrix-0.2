#!/bin/bash

xdg-icon-resource install --context apps --size 64 icons/mrtrix.png application-x-mrtrix
xdg-icon-resource install --context apps --size 64 icons/mrtrix.png application-x-mrtrix-gz
xdg-icon-resource install --context apps --size 64 icons/mrtrix.png application-x-nifti
xdg-icon-resource install --context apps --size 64 icons/mrtrix.png application-x-nifti-gz
xdg-icon-resource install --context apps --size 64 icons/mrtrix.png application-x-mgh
xdg-icon-resource install --context apps --size 64 icons/mrtrix.png application-x-mgz
xdg-icon-resource install --context apps --size 64 icons/mrtrix.png application-x-analyze
xdg-icon-resource install --context apps --size 64 icons/tracks.png application-x-mrtrix-tracks

xdg-mime install mrtrix-mime.xml
xdg-desktop-menu install mrtrix-mrview.desktop

