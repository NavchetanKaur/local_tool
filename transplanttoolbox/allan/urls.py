"""transplanttoolbox URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.conf.urls import url
from django.conf import settings
from allan import views, views_a, views_al, views_gl, views_ac, views_reverse
from rest_framework_swagger.views import get_swagger_view

schema_view = get_swagger_view(title='Web services interface information for Transplanttoolbox')

from django.urls import path

urlpatterns = [
    path('', views.home, name="home"),
    path('home_1', views.home_1, name="home_1"),
    path('cpra_message', views.cpra_maintenance_message, name="cpra_message"),
    path('license', views.license, name="license"),
    path('allele', views.allele, name="allele"),
    path('lists', views.allele_list, name="lists"),
    path('gl_string', views.gl_string, name="gl_string"),
    path('codes', views.allele_codes, name="codes"),
    path('convert/', views.convert, name="convert"),
   	path('convert_2/', views.convert_2, name="convert_2"),
    path('convert_3/', views.convert_3, name="convert_3"),
    path('convert_4/', views.convert_4, name="convert_4"),
    path('reverse/', views.reverse, name="reverse"),
    path('convert_5/', views.convert_5, name="convert_5"),
    path('single_allele/', views_a.AlleleApiView.as_view(), name='asaw'),
    path('array/', views_al.AlleleListApiView.as_view()),
    path('gls/', views_gl.GLstringApiView.as_view()),
    path('macs/', views_ac.AlleleCodesApiView.as_view()),
    path('reverse_mapping/', views_reverse.AlleleMappingApiView.as_view()),
    path('services', schema_view, name="allan_services"),
]
