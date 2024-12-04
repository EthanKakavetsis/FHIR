import React, { MouseEventHandler, useCallback, useState, useRef } from "react";
import axios from "axios"
import SortableTable from "./SortableTable";

const baseURLFSV = 'https://fhir-gen-ops.herokuapp.com/subject-operations/genotype-operations/$find-subject-variants?subject='
const baseURLFC = 'https://fhir-gen-ops.herokuapp.com/utilities/get-feature-coordinates?gene='

const config = {
    headers: {
        "Access-Control-Allow-Origin": "*",
        "Access-Control-Allow-Methods": "GET,PUT,POST,DELETE,PATCH,OPTIONS"
    }
}

type VariantRow = {
    spdi: string,
    dnaChangeType: string,
    sourceClass: string,
    allelicState: string,
    molecImpact: string,
    alleleFreq: number,
}

type FSVResponse = {
    parameter: Array<{
        name: string,
        part: Array<FSVParameter>
    }>,
    resourceType: string
}

type FHIRCoding = [{
    code: string,
    system: string,
    display?: string
}]

type FSVVariantResource = {
    category: Array<FHIRCoding>,
    code: FHIRCoding,
    component: Array<{
        code: { coding: FHIRCoding },
        valueCodeableConcept?: { coding: FHIRCoding },
        valueString?: string,
        interpretation?: { text: string },
        valueQuantity?: { code: string, system: string, value: number },
        valueRange?: { low: { value: number } }
    }>,
    id: string,
    meta: { profile: Array<string> },
    resourceType: string,
    status: string,
    subject: { reference: string },
    valueCodeableConcept: FHIRCoding,
    text?: string,
}

type FSVParameter = {
    name: string,
    resource?: FSVVariantResource,
    valueString?: string,
    valueBoolean?: boolean,
}

type variant = {
    id: string,
}

function translateComponents(resource: FSVVariantResource) {
    let variantRow: VariantRow = {
        spdi: "",
        dnaChangeType: "",
        sourceClass: "",
        allelicState: "",
        molecImpact: "",
        alleleFreq: NaN,
    }
    resource.component.map((component) => {
        if (component.code.coding[0].code == "48002-0") {
            if (component.valueCodeableConcept && component.valueCodeableConcept.coding[0].display) {
                variantRow.sourceClass = component.valueCodeableConcept.coding[0].display
            }
        } else if (component.code.coding[0].code == "53034-5") {
            if (component.valueCodeableConcept && component.valueCodeableConcept.coding[0].display) {
                variantRow.allelicState = component.valueCodeableConcept.coding[0].display
            }
        } else if (component.code.coding[0].code == "81252-9") {
            if (component.valueCodeableConcept && component.valueCodeableConcept.coding[0].display) {
                variantRow.spdi = component.valueCodeableConcept.coding[0].display
            }
        } else if (component.code.coding[0].code == "81258-6") {
            if (component.valueQuantity && component.valueQuantity.value) {
                variantRow.alleleFreq = component.valueQuantity.value
            }
        } else if (component.code.coding[0].code == "81258-6") {
            if (component.valueQuantity && component.valueQuantity.value) {
                variantRow.alleleFreq = component.valueQuantity.value
            }
        } else if (component.code.coding[0].code == "molecular-consequence") {
            if (component.interpretation) {
                variantRow.molecImpact = component.interpretation.text
            }
        }
    })
    return variantRow
}

function translateFHIRResponse(response: FSVResponse) {
    let geneTable: Array<VariantRow> = []
    let paramArray = response.parameter[0].part;
    paramArray.map((param) => {
        //If the part is not a variant resource, return
        if (param.name != 'variant' || param.resource == null) {
            return
        }

        let variantRow: VariantRow = translateComponents(param.resource)
        geneTable.push(variantRow)
    })

    return geneTable
}

function GeneButton({ patientID, gene, addAnnFlag, callback }:
    {
        patientID: string,
        gene: string,
        addAnnFlag: boolean,
        callback: (geneData: Array<VariantRow>) => void
    }) {

    // const [APIStatus, setAPIStatus] = useState("red")
    var APIStatus = "red"
    const range = useRef("")
    var geneData: Array<VariantRow>

    console.log("In gene handler")
    console.log("Gene:" + gene)

    console.log("PatientID: " + patientID)

    if (!patientID) {
        return null
    }

    async function getFC() {
        if (range.current != "") {
            return
        }
        type FeatureCoordinates = [{
            MANE: string[],
            build37Coordinates: string,
            build38Coordinates: string,
            geneId: string,
            geneLink: string,
            geneSymbol: string,
            transcripts: string[]
        }]

        let url = baseURLFC + gene
        // urlAppender = urlAppender.replaceAll("/", "@")
        // urlAppender = urlAppender.replaceAll("$", "@!abab@!")
        // urlAppender = urlAppender.replaceAll("?", "!")


        // let url = `http://127.0.0.1:5000/${urlAppender}`

        let response: any

        response = await fetch(url)
        var responseJson = await response.json() as FeatureCoordinates
        if (responseJson instanceof Error) {
            console.log('It is an error!');
        }
        else {
            console.log(responseJson);
            // setAPIStatus("yellow")
            APIStatus = 'yellow'
            // setArticleDict(JSON.parse(responseJson));
        }

        range.current = responseJson[0].build38Coordinates
    }

    async function getFSV() {
        if (range.current == "") {
            return
        }

        let url = baseURLFSV + patientID + '&ranges=' + range.current + '&includeVariants=true'
        // urlAppender = urlAppender.replaceAll("/", "@")
        // urlAppender = urlAppender.replaceAll("$", "@!abab@!")
        // urlAppender = urlAppender.replaceAll("?", "!")

        // let url = `http://127.0.0.1:5000/${urlAppender}`
        console.log(url)

        let response = await fetch(url)
        var fsvResponseJson = await response.json() as FSVResponse
        if (fsvResponseJson instanceof Error) {
            console.log('It is an error!');
        }
        else {
            console.log(fsvResponseJson.parameter[0]);
            geneData = translateFHIRResponse(fsvResponseJson);
            // setAPIStatus("green")
            APIStatus = 'green'
        }
    }

    async function getFHIRResponse() {
        await getFC()
        await getFSV()
    }

    getFHIRResponse()

    return (
        <div>
            <button style={{ color: APIStatus }} onClick={() => {
                if (APIStatus == 'green') {
                    console.log("Pushed Gene Button!")
                    callback(geneData)
                }
            }}>{gene}</button>

        </div >
    )
}
export default GeneButton
